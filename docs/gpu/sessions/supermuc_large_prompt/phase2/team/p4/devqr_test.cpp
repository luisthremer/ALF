// Standalone correctness + timing test of the device QR (Libraries/Modules/alf_la_sycl_devqr.hpp) against host MKL ZGEQRF.
// Build: icpx -fsycl -O2 -std=c++17 -I<repo>/Libraries/Modules -qmkl=parallel devqr_test.cpp -o devqr_test
// Run:   SYCL_UR_USE_LEVEL_ZERO_V2=1 ZE_AFFINITY_MASK=0 ./devqr_test [nb] [timing reps]      env DEVQR_PLANT=1: planted defect tau(0)*=1+1e-12
#include <complex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <sycl/sycl.hpp>
#include <oneapi/mkl/lapack.hpp>
#include "alf_la_sycl_devqr.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <random>
#include <vector>
using cd = std::complex<double>;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static void gen(int type, int m, int n, std::vector<cd> &A, unsigned seed) {
  std::mt19937_64 g(seed);
  std::normal_distribution<double> u(0, 1);
  A.assign((size_t)m * n, 0);
  for (auto &x : A) x = cd(u(g), u(g));
  if (type == 1 || type == 2) {   // column scaling 1e-12..1e12: sorted decreasing (after QR pre-pivot) / random order
    std::vector<int> p(n);
    for (int j = 0; j < n; ++j) p[j] = j;
    if (type == 2) std::shuffle(p.begin(), p.end(), g);
    for (int j = 0; j < n; ++j) {
      double s = std::pow(10.0, 12.0 - 24.0 * p[j] / std::max(1, n - 1));
      for (int i = 0; i < m; ++i) A[i + (size_t)j * m] *= s;
    }
  } else if (type == 3) {   // genuinely ill-conditioned: A = X * diag(10^(12..-12)) * Y (cond ~1e24), orthogonal-ish X, Y mixing
    std::vector<cd> X((size_t)m * n), Y((size_t)n * n), T((size_t)m * n);
    for (auto &x : X) x = cd(u(g), u(g));
    for (auto &x : Y) x = cd(u(g), u(g));
    for (int j = 0; j < n; ++j) { double s = std::pow(10.0, 12.0 - 24.0 * j / std::max(1, n - 1)); for (int i = 0; i < m; ++i) X[i + (size_t)j * m] *= s; }
    cd one(1, 0), zero(0, 0);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, n, &one, X.data(), m, Y.data(), n, &zero, A.data(), m);
  } else if (type == 4) {   // extreme scales 1e-150..1e150 (overflow/underflow handling) plus one 1e-300 column
    for (int j = 0; j < n; ++j) {
      double s = std::pow(10.0, 150.0 - 300.0 * j / std::max(1, n - 1));
      for (int i = 0; i < m; ++i) A[i + (size_t)j * m] *= s;
    }
  }
}
// returns max |Q^H Q - I|
static double orth(int m, int n, const std::vector<cd> &F, const std::vector<cd> &tau) {
  std::vector<cd> Q(F);
  std::vector<cd> t(tau);
  MKL_INT k = std::min(m, n), info, mm = m;
  MKL_INT ncol = m >= n ? n : m;   // Q is m x ncol
  std::vector<cd> Qm((size_t)m * ncol);
  for (int j = 0; j < ncol; ++j) for (int i = 0; i < m; ++i) Qm[i + (size_t)j * m] = F[i + (size_t)j * m];
  cd wq; MKL_INT lw = -1; zungqr(&mm, &ncol, &k, Qm.data(), &mm, t.data(), &wq, &lw, &info); lw = (MKL_INT)wq.real();
  std::vector<cd> w(lw); zungqr(&mm, &ncol, &k, Qm.data(), &mm, t.data(), w.data(), &lw, &info);
  std::vector<cd> G((size_t)ncol * ncol);
  cd one(1, 0), zero(0, 0);
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, ncol, ncol, m, &one, Qm.data(), m, Qm.data(), m, &zero, G.data(), ncol);
  double mx = 0;
  for (int j = 0; j < ncol; ++j) for (int i = 0; i < ncol; ++i) mx = std::max(mx, std::abs(G[i + (size_t)j * ncol] - (i == j ? 1.0 : 0.0)));
  return mx;
}

// copies through compute kernels and pinned host USM (no copy engine): DEVQR_KCOPY=1
// NOTE (SuperMUC-NG Phase 2, SYCL_UR_USE_LEVEL_ZERO_V2=1, which standalone programs need): a device->host copy is complete only some
// ms AFTER queue.wait()/event.wait() return (sanity.cpp: ~1% of a 1 MB readback still zero when checked at once, 0 after a 20 ms sleep),
// so every readback below sleeps 30 ms. Kernel completion itself is not affected.
static void kcopy(sycl::queue &q, void *dst, const void *src, size_t bytes, bool dst_is_dev) {
  static const bool k = std::getenv("DEVQR_KCOPY") != nullptr;
  if (!k) { q.memcpy(dst, src, bytes).wait(); if (!dst_is_dev) std::this_thread::sleep_for(std::chrono::milliseconds(30)); return; }
  size_t nd = bytes / 8;
  void *stage = sycl::malloc_host(bytes, q);
  if (dst_is_dev) { std::memcpy(stage, src, bytes); double *d = (double *)dst; const double *s = (const double *)stage; q.parallel_for(nd, [=](sycl::id<1> i) { d[i] = s[i]; }).wait(); }
  else { double *d = (double *)stage; const double *s = (const double *)src; q.parallel_for(nd, [=](sycl::id<1> i) { d[i] = s[i]; }).wait(); std::memcpy(dst, stage, bytes); }
  sycl::free(stage, q);
}
int main(int argc, char **argv) {
  int nb = argc > 1 ? atoi(argv[1]) : 16;
  int reps = argc > 2 ? atoi(argv[2]) : 20;
  bool plant = std::getenv("DEVQR_PLANT") != nullptr;
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order()};
  std::printf("device %s nb=%d plant=%d mkl_threads=%d\n", q.get_device().get_info<sycl::info::device::name>().c_str(), nb, (int)plant, mkl_get_max_threads());
  struct Case { int m, n, type; };
  std::vector<Case> cases;
  for (int type = 0; type <= 4; ++type) for (int n : {36, 144, 576, 1024}) cases.push_back({n, n, type});
  cases.push_back({576, 144, 0}); cases.push_back({144, 576, 0}); cases.push_back({1000, 37, 1}); cases.push_back({37, 1000, 1});
  cases.push_back({5, 5, 0}); cases.push_back({1, 1, 0}); cases.push_back({17, 33, 2}); cases.push_back({577, 577, 3});
  if (argc > 3) cases.clear();   // timing only
  if (const char *ce = std::getenv("DEVQR_CASES")) {   // e.g. "576x16:0,576x32:0"; each repeated 4x
    cases.clear();
    for (const char *p = ce; *p;) { int m, n, t; int c; if (sscanf(p, "%dx%d:%d%n", &m, &n, &t, &c) != 3) break; for (int r = 0; r < 4; ++r) cases.push_back({m, n, t}); p += c; if (*p == ',') ++p; }
  }
  int nfail = 0;
  for (auto &cs : cases) {
    int m = cs.m, n = cs.n, k = std::min(m, n);
    std::vector<cd> A0; gen(cs.type, m, n, A0, 12345 + m + n);
    // host reference
    std::vector<cd> H(A0), th(k);
    { MKL_INT info, mm = m, nn = n, lw = -1; cd wq; zgeqrf(&mm, &nn, H.data(), &mm, th.data(), &wq, &lw, &info); lw = (MKL_INT)wq.real();
      std::vector<cd> w(lw); zgeqrf(&mm, &nn, H.data(), &mm, th.data(), w.data(), &lw, &info); }
    // device
    int64_t ws = alfdevqr::scratch_size(m, n, nb);
    cd *dA = sycl::malloc_device<cd>((size_t)m * n, q), *dt = sycl::malloc_device<cd>(k, q), *dw = sycl::malloc_device<cd>(ws, q);
    kcopy(q, dA, A0.data(), sizeof(cd) * (size_t)m * n, true);
    sycl::event ev;
    int rc = 0;
    if (std::getenv("DEVQR_USEMKL")) { int64_t w2 = oneapi::mkl::lapack::geqrf_scratchpad_size<cd>(q, m, n, m); cd *dw2 = sycl::malloc_device<cd>(w2, q); oneapi::mkl::lapack::geqrf(q, m, n, dA, m, dt, dw2, w2).wait_and_throw(); sycl::free(dw2, q); }
    else rc = alfdevqr::geqrf(q, m, n, dA, m, dt, dw, ws, nb, &ev);
    q.wait_and_throw();
    std::vector<cd> D((size_t)m * n), td(k);
    kcopy(q, D.data(), dA, sizeof(cd) * (size_t)m * n, false); kcopy(q, td.data(), dt, sizeof(cd) * k, false);
    if (plant) td[0] *= (1 + 1e-12);
    // host MKL ZGEQRF leaves H(k)=I (tau=0, R(k,k) complex) when the last reflector acts on ONE row (m<=n); reference LAPACK and the
    // device oneMKL ZGEQRF (what ALF uses today) make it real. Such a row k-1 of R / tau(k-1) is excluded from the comparison with the host.
    const bool skipLast = (m <= n) && std::abs(th[k - 1]) == 0 && std::abs(td[k - 1]) != 0;
    // deviations, column-wise normalised
    double dR = 0, dV = 0, dT = 0;
    for (int j = 0; j < n; ++j) {
      double nr = 0, nv = 0, er = 0, ev2 = 0;
      for (int i = 0; i < m; ++i) {
        cd h = H[i + (size_t)j * m], d = D[i + (size_t)j * m];
        if (skipLast && i == k - 1) continue;
        if (i <= j) { nr = std::max(nr, std::abs(h)); er = std::max(er, std::abs(h - d)); }
        else { nv = std::max(nv, std::abs(h)); ev2 = std::max(ev2, std::abs(h - d)); }
      }
      // R normalised by the column's largest R entry; reflector by its largest component
      if (nr > 0) { if (std::getenv("DEVQR_DBG") && er / nr > 1e-13) { for (int i = 0; i <= std::min(j, m - 1); ++i) if (std::abs(H[i + (size_t)j * m] - D[i + (size_t)j * m]) > 1e-13 * nr) { std::printf("  bad R(%d,%d) host %g,%g dev %g,%g\n", i, j, H[i + (size_t)j * m].real(), H[i + (size_t)j * m].imag(), D[i + (size_t)j * m].real(), D[i + (size_t)j * m].imag()); break; } } dR = std::max(dR, er / nr); }
      if (nv > 0) { if (std::getenv("DEVQR_DBG") && ev2 / nv > 1e-13) { std::printf("  bad V col %d nv=%g err=%g rows:", j, nv, ev2); int cnt = 0; for (int i = j + 1; i < m && cnt < 40; ++i) if (std::abs(H[i + (size_t)j * m] - D[i + (size_t)j * m]) > 1e-13 * nv) { std::printf(" %d(%.3g)", i, std::abs(D[i + (size_t)j * m])); ++cnt; } std::printf("\n"); } dV = std::max(dV, ev2 / nv); }
    }
    if (std::getenv("DEVQR_DBG") && (dV > 1e-6 || dR > 1e-6) && n <= 200 && !std::getenv("DEVQR_NOMAP")) {   // mismatch map in 16x16 blocks: '.' ok, 'Z' dev block all zero, 'x' differs
      for (int bi = 0; bi < (m + 15) / 16; ++bi) { for (int bj = 0; bj < (n + 15) / 16; ++bj) { int nz = 0, nd = 0, nt = 0; for (int i = bi * 16; i < std::min(m, bi * 16 + 16); ++i) for (int j = bj * 16; j < std::min(n, bj * 16 + 16); ++j) { ++nt; if (std::abs(D[i + (size_t)j * m]) == 0) ++nz; if (std::abs(D[i + (size_t)j * m] - H[i + (size_t)j * m]) > 1e-9 * (1 + std::abs(H[i + (size_t)j * m]))) ++nd; } std::putchar(nd == 0 ? '.' : (nz == nt ? 'Z' : 'x')); } std::putchar('\n'); }
    }
    for (int j = 0; j < k - (skipLast ? 1 : 0); ++j) if (std::abs(th[j]) > 0) dT = std::max(dT, std::abs(th[j] - td[j]) / std::abs(th[j]));
    // full comparison with the device oneMKL ZGEQRF (the routine devqr replaces)
    double dRm = 0, dVm = 0, dTm = 0;
    if (std::getenv("DEVQR_VSMKL")) {   // NOT reliable under SYCL_UR_USE_LEVEL_ZERO_V2=1 (hybrid oneMKL ZGEQRF returns wrong results there): off by default
      int64_t w2 = oneapi::mkl::lapack::geqrf_scratchpad_size<cd>(q, m, n, m); cd *dw2 = sycl::malloc_device<cd>(w2, q), *dB = sycl::malloc_device<cd>((size_t)m * n, q), *dt2 = sycl::malloc_device<cd>(k, q);
      kcopy(q, dB, A0.data(), sizeof(cd) * (size_t)m * n, true);
      oneapi::mkl::lapack::geqrf(q, m, n, dB, m, dt2, dw2, w2).wait_and_throw();
      std::vector<cd> Bm((size_t)m * n), tm(k); kcopy(q, Bm.data(), dB, sizeof(cd) * (size_t)m * n, false); kcopy(q, tm.data(), dt2, sizeof(cd) * k, false);
      for (int j = 0; j < n; ++j) { double nr = 0, nv = 0, er = 0, ev2 = 0; for (int i = 0; i < m; ++i) { cd h = Bm[i + (size_t)j * m], d = D[i + (size_t)j * m]; if (i <= j) { nr = std::max(nr, std::abs(h)); er = std::max(er, std::abs(h - d)); } else { nv = std::max(nv, std::abs(h)); ev2 = std::max(ev2, std::abs(h - d)); } } if (nr > 0) dRm = std::max(dRm, er / nr); if (nv > 0) dVm = std::max(dVm, ev2 / nv); }
      for (int j = 0; j < k; ++j) if (std::abs(tm[j]) > 0) dTm = std::max(dTm, std::abs(tm[j] - td[j]) / std::abs(tm[j]));
      sycl::free(dw2, q); sycl::free(dB, q); sycl::free(dt2, q);
    }
    double oh = 0, od = 0;
    if (cs.type != 4) { oh = orth(m, n, H, th); od = orth(m, n, D, td); }   // type 4: columns beyond double range of Q^H Q? (Q unaffected; skipped only for speed)
    else { oh = orth(m, n, H, th); od = orth(m, n, D, td); }
    // type 3 (cond ~1e24 > 1/eps): the trailing reflectors are rounding noise, any two correct QRs differ there (R and Q stay equivalent):
    // informational for R/V/tau, gated only on the orthogonality of Q
    const bool strict = cs.type != 3;
    bool ok = rc == 0 && (!strict || (dR <= 1e-13 && dV <= 1e-13 && dT <= 1e-13)) && od <= 2 * oh + 1e-15;
    if (!ok) ++nfail;
    std::printf("%s m=%4d n=%4d type=%d rc=%d  dR=%.2e dV=%.2e dtau=%.2e  orth host=%.2e dev=%.2e%s %s\n", ok ? (strict ? "PASS" : "PASS(info)") : "FAIL", m, n, cs.type, rc, dR, dV, dT, oh, od, skipLast ? " [last row skipped]" : "", std::getenv("DEVQR_VSMKL") ? "(vs oneMKL-dev: see DEVQR_VSMKL)" : "");
    sycl::free(dA, q); sycl::free(dt, q); sycl::free(dw, q);
  }
  std::printf("correctness: %d failed\n", nfail);
  if (std::getenv("DEVQR_CASES")) return nfail ? 1 : 0;
  if (argc > 3) std::printf("(timing only)\n");
  // timing: device QR vs oneMKL (hybrid) per call, square random
  for (int n : {576, 1024}) {
    int tty = std::getenv("DEVQR_TTYPE") ? atoi(std::getenv("DEVQR_TTYPE")) : 0;
    std::vector<cd> A0; gen(tty, n, n, A0, 7);
    int64_t ws = std::max<int64_t>(alfdevqr::scratch_size(n, n, nb), oneapi::mkl::lapack::geqrf_scratchpad_size<cd>(q, n, n, n));
    cd *A = sycl::malloc_device<cd>((size_t)n * n, q), *A00 = sycl::malloc_device<cd>((size_t)n * n, q), *t = sycl::malloc_device<cd>(n, q), *w = sycl::malloc_device<cd>(ws, q);
    q.memcpy(A00, A0.data(), sizeof(cd) * (size_t)n * n).wait();
    for (int mode = 0; mode < 2; ++mode) {
      double tt = 0, tmin = 1e9;
      for (int r = 0; r < reps + 2; ++r) {
        q.memcpy(A, A00, sizeof(cd) * (size_t)n * n).wait();
        double t0 = now();
        if (mode == 0) { sycl::event e; alfdevqr::geqrf(q, n, n, A, n, t, w, ws, nb, &e); e.wait_and_throw(); }
        else oneapi::mkl::lapack::geqrf(q, n, n, A, n, t, w, ws).wait_and_throw();
        double dt = now() - t0;
        if (r >= 2) { tt += dt; tmin = std::min(tmin, dt); }
      }
      std::printf("TIMING(V2: submission time only, NOT kernel time) type=%d N=%d %-8s mean %8.3f ms  min %8.3f ms\n", tty, n, mode == 0 ? "devqr" : "oneMKL", 1e3 * tt / reps, 1e3 * tmin);
    }
  }
  return nfail ? 1 : 0;
}
