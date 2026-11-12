// Standalone test of alfhq_geqrf / alfhq_unmqr against hipSOLVER zgeqrf / zunmqr.
// build: hipcc -O2 --offload-arch=gfx90a -I$ROCM_INC qr_test.cpp ../../../Libraries/Modules/alf_hip_qr.hip.cpp -lhipsolver
#include <hip/hip_runtime.h>
#include <hipsolver/hipsolver.h>
#include <complex>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
typedef std::complex<double> zc;
extern "C" int alfhq_geqrf(int, int, void*, int, void*);
extern "C" int alfhq_unmqr(int, int, int, int, int, void*, int, void*, void*, int);
static double rnd() { return rand() / (double)RAND_MAX - 0.5; }
static double nrm(const std::vector<zc>& a) { double s = 0; for (auto& x : a) s += std::norm(x); return sqrt(s); }
static double diffn(const std::vector<zc>& a, const std::vector<zc>& b) { double s = 0; for (size_t i = 0; i < a.size(); i++) s += std::norm(a[i] - b[i]); return sqrt(s); }
static hipsolverHandle_t h;
static zc *dW = nullptr; static int *dinfo; static int Lw = 0;
static void need(int l) { if (l > Lw) { if (dW) hipFree(dW); hipMalloc(&dW, (size_t)l * 16); Lw = l; } }
int main(int argc, char** argv) {
  hipsolverCreate(&h); hipMalloc(&dinfo, 4);
  struct Cs { int m, n; int graded; } cases[] = {{576, 576, 0}, {576, 576, 1}, {300, 200, 0}, {200, 300, 0}, {1000, 1000, 1}, {64, 64, 0}, {17, 17, 0}, {1, 1, 0}, {5, 3, 0}, {1024, 1024, 0}, {2304, 2304, 1}, {1500, 700, 0}};
  int fails = 0;
  for (auto cs : cases) {
    int m = cs.m, n = cs.n, k = std::min(m, n);
    std::vector<zc> A((size_t)m * n);
    for (int j = 0; j < n; j++) { double sc = cs.graded ? pow(10.0, -8.0 * j / std::max(1, n - 1)) : 1.0; for (int i = 0; i < m; i++) A[i + (size_t)j * m] = zc(rnd(), rnd()) * sc; }
    if (cs.m == 5) A[0] = zc(0, 0);                     // exercise the tau = 0 / zero-column paths
    zc *d1, *d2, *t1, *t2; size_t sz = (size_t)m * n * 16;
    hipMalloc(&d1, sz); hipMalloc(&d2, sz); hipMalloc(&t1, k * 16); hipMalloc(&t2, k * 16);
    hipMemcpy(d1, A.data(), sz, hipMemcpyHostToDevice); hipMemcpy(d2, A.data(), sz, hipMemcpyHostToDevice);
    int lw; hipsolverZgeqrf_bufferSize(h, m, n, (hipDoubleComplex*)d1, m, &lw); need(lw);
    hipsolverZgeqrf(h, m, n, (hipDoubleComplex*)d1, m, (hipDoubleComplex*)t1, (hipDoubleComplex*)dW, lw, dinfo); hipDeviceSynchronize();
    int rc = alfhq_geqrf(m, n, d2, m, t2);
    std::vector<zc> A1((size_t)m * n), A2((size_t)m * n), T1(k), T2(k);
    hipMemcpy(A1.data(), d1, sz, hipMemcpyDeviceToHost); hipMemcpy(A2.data(), d2, sz, hipMemcpyDeviceToHost);
    hipMemcpy(T1.data(), t1, k * 16, hipMemcpyDeviceToHost); hipMemcpy(T2.data(), t2, k * 16, hipMemcpyDeviceToHost);
    // column-wise relative difference of the factors (columns can differ in scale for graded matrices)
    double worst = 0; for (int j = 0; j < n; j++) { double s = 0, d = 0; for (int i = 0; i < m; i++) { s += std::norm(A1[i + (size_t)j * m]); d += std::norm(A1[i + (size_t)j * m] - A2[i + (size_t)j * m]); } worst = std::max(worst, sqrt(d) / std::max(sqrt(s), 1e-300)); }
    double wt = 0; for (int j = 0; j < k; j++) wt = std::max(wt, std::abs(T1[j] - T2[j]) / std::max(std::abs(T1[j]), 1e-300));
    // residual: Q^H A (from alfhq_unmqr, left C) must equal R (upper triangle of the factors) -> check || Q^H A - R || / ||A||
    hipMemcpy(d1, A.data(), sz, hipMemcpyHostToDevice);
    rc |= alfhq_unmqr(0, 1, m, n, k, d2, m, t2, d1, m);
    std::vector<zc> QA((size_t)m * n); hipMemcpy(QA.data(), d1, sz, hipMemcpyDeviceToHost);
    double res = 0, an = nrm(A); for (int j = 0; j < n; j++) for (int i = 0; i < m; i++) { zc r = (i <= j && i < k) ? A2[i + (size_t)j * m] : zc(0, 0); res += std::norm(QA[i + (size_t)j * m] - r); }
    // all four unmqr variants vs hipSOLVER on a random C
    double wu = 0; const char* nm[4] = {"LN", "LC", "RN", "RC"};
    double wv[4];
    for (int v = 0; v < 4; v++) {
      int side = v / 2 == 0 ? 0 : 1; int tr = (v % 2 == 0 && side == 0) ? 0 : (side == 0 ? 2 : (v == 2 ? 0 : 2));
      if (side == 0 && v == 1) tr = 2; if (side == 1) tr = (v == 2) ? 0 : 2;
      int cm = m, cn = n; if (side == 1) { cm = 3 * 7 + 5; cn = m; if (cn < k) continue; }
      if (side == 0) { cm = m; cn = 37; }
      std::vector<zc> C((size_t)cm * cn); for (auto& x : C) x = zc(rnd(), rnd());
      zc *c1, *c2; hipMalloc(&c1, C.size() * 16); hipMalloc(&c2, C.size() * 16);
      hipMemcpy(c1, C.data(), C.size() * 16, hipMemcpyHostToDevice); hipMemcpy(c2, C.data(), C.size() * 16, hipMemcpyHostToDevice);
      // hipSOLVER reference uses its own factors d1' (recompute into a fresh copy so tau/A are consistent)
      zc *d3, *t3; hipMalloc(&d3, sz); hipMalloc(&t3, k * 16); hipMemcpy(d3, A.data(), sz, hipMemcpyHostToDevice);
      hipsolverZgeqrf(h, m, n, (hipDoubleComplex*)d3, m, (hipDoubleComplex*)t3, (hipDoubleComplex*)dW, lw, dinfo);
      hipsolverSideMode_t sd = side == 0 ? HIPSOLVER_SIDE_LEFT : HIPSOLVER_SIDE_RIGHT;
      hipsolverOperation_t op = tr == 0 ? HIPSOLVER_OP_N : HIPSOLVER_OP_C;
      int l2; hipsolverZunmqr_bufferSize(h, sd, op, cm, cn, k, (hipDoubleComplex*)d3, m, (hipDoubleComplex*)t3, (hipDoubleComplex*)c1, cm, &l2); need(l2);
      hipsolverZunmqr(h, sd, op, cm, cn, k, (hipDoubleComplex*)d3, m, (hipDoubleComplex*)t3, (hipDoubleComplex*)c1, cm, (hipDoubleComplex*)dW, l2, dinfo); hipDeviceSynchronize();
      rc |= alfhq_unmqr(side, tr, cm, cn, k, d3, m, t3, c2, cm);
      std::vector<zc> r1(C.size()), r2(C.size()); hipMemcpy(r1.data(), c1, C.size() * 16, hipMemcpyDeviceToHost); hipMemcpy(r2.data(), c2, C.size() * 16, hipMemcpyDeviceToHost);
      wv[v] = diffn(r1, r2) / nrm(r1); wu = std::max(wu, wv[v]);
      hipFree(c1); hipFree(c2); hipFree(d3); hipFree(t3);
    }
    bool ok = rc == 0 && worst < 1e-10 && res / an < 1e-13 && wu < 1e-13;
    printf("%s m=%4d n=%4d graded=%d rc=%d  factors-vs-hipsolver(max col rel) %.1e  tau %.1e  ||Q^H A-R||/||A|| %.1e  unmqr(LN,LC,RN,RC) rel %.1e %.1e %.1e %.1e\n",
           ok ? "PASS" : "FAIL", m, n, cs.graded, rc, worst, wt, res / an, wv[0], wv[1], wv[2], wv[3]);
    if (!ok) fails++;
    hipFree(d1); hipFree(d2); hipFree(t1); hipFree(t2);
  }
  // timing
  for (int n : {576, 1024, 2304}) {
    size_t sz = (size_t)n * n * 16; std::vector<zc> A((size_t)n * n); for (auto& x : A) x = zc(rnd(), rnd());
    zc *d0, *d1, *d2, *t1, *t2, *c; hipMalloc(&d0, sz); hipMalloc(&d1, sz); hipMalloc(&d2, sz); hipMalloc(&c, sz); hipMalloc(&t1, n * 16); hipMalloc(&t2, n * 16);
    hipMemcpy(d0, A.data(), sz, hipMemcpyHostToDevice);
    int lw, l2; hipsolverZgeqrf_bufferSize(h, n, n, (hipDoubleComplex*)d0, n, &lw); hipsolverZunmqr_bufferSize(h, HIPSOLVER_SIDE_LEFT, HIPSOLVER_OP_N, n, n, n, (hipDoubleComplex*)d0, n, (hipDoubleComplex*)t1, (hipDoubleComplex*)c, n, &l2); need(std::max(lw, l2));
    auto tm = [&](auto f) { double t = 0; int reps = 10; for (int r = 0; r <= reps; r++) { hipMemcpy(d1, d0, sz, hipMemcpyDeviceToDevice); hipMemcpy(c, d0, sz, hipMemcpyDeviceToDevice); hipDeviceSynchronize(); auto t0 = std::chrono::steady_clock::now(); f(); hipDeviceSynchronize(); if (r) t += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); } return 1e3 * t / reps; };
    double g1 = tm([&] { hipsolverZgeqrf(h, n, n, (hipDoubleComplex*)d1, n, (hipDoubleComplex*)t1, (hipDoubleComplex*)dW, lw, dinfo); });
    double g2 = tm([&] { alfhq_geqrf(n, n, d1, n, t2); });
    hipMemcpy(d2, d0, sz, hipMemcpyDeviceToDevice); alfhq_geqrf(n, n, d2, n, t2);
    double u1 = tm([&] { hipsolverZunmqr(h, HIPSOLVER_SIDE_LEFT, HIPSOLVER_OP_N, n, n, n, (hipDoubleComplex*)d2, n, (hipDoubleComplex*)t2, (hipDoubleComplex*)c, n, (hipDoubleComplex*)dW, l2, dinfo); });
    double u2 = tm([&] { alfhq_unmqr(0, 0, n, n, n, d2, n, t2, c, n); });
    double u3 = tm([&] { alfhq_unmqr(1, 0, n, n, n, d2, n, t2, c, n); });
    printf("n=%d  geqrf: hipsolver %.2f ms, alf %.2f ms | unmqr LN: hipsolver %.2f ms, alf %.2f ms (RN alf %.2f ms)\n", n, g1, g2, u1, u2, u3);
  }
  printf("%s\n", fails ? "SOME FAILED" : "ALL PASSED");
  return fails != 0;
}
