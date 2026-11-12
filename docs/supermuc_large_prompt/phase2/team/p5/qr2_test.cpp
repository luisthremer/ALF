// Phase 5 standalone test: T-free device QR (geqrf2) and device ZUNMQR (unmqr2) of Libraries/Modules/alf_la_sycl_qr2.hpp against host MKL.
// Build: icpx -fsycl -O2 -std=c++17 -I<repo>/Libraries/Modules -qmkl=parallel qr2_test.cpp -o qr2_test      [AOT: -fsycl-targets=spir64_gen -Xs "-device pvc"]
// Run:   SYCL_UR_USE_LEVEL_ZERO_V2=1 ZE_AFFINITY_MASK=0 ./qr2_test <mode> [nb]
//   mode "correct": Gate A (geqrf: R, reflectors, tau, orthogonality; unmqr: L/R x N/C vs host zunmqr, Q^H Q = I). DEVQR_PLANT=1|2|3 plants a defect that must FAIL:
//        1 tau(0) *= 1+1e-12 (geqrf output), 2 wrong conj (device unmqr asked for the other trans), 3 reflector scaled by 1+1e-12 (unmqr input)
//   mode "time"   : true device time per call (chain of calls closed by a kernel that sets a host-USM flag, host spins; NOT queue.wait(): see p4 STATUS)
//   mode "panel"  : time of one panel kernel vs M and jb
#include <complex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <sycl/sycl.hpp>
#include "alf_la_sycl_qr2.hpp"
#include <chrono>
#include <thread>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <string>
using cd = std::complex<double>;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static void gen(int type, int m, int n, std::vector<cd> &A, unsigned seed) {
  std::mt19937_64 g(seed);
  std::normal_distribution<double> u(0, 1);
  A.assign((size_t)m * n, 0);
  for (auto &x : A) x = cd(u(g), u(g));
  if (type == 1 || type == 2) {   // column scaling 1e-12..1e12: sorted decreasing / random order
    std::vector<int> p(n);
    for (int j = 0; j < n; ++j) p[j] = j;
    if (type == 2) std::shuffle(p.begin(), p.end(), g);
    for (int j = 0; j < n; ++j) { double s = std::pow(10.0, 12.0 - 24.0 * p[j] / std::max(1, n - 1)); for (int i = 0; i < m; ++i) A[i + (size_t)j * m] *= s; }
  } else if (type == 4) {   // extreme scales 1e-150..1e150
    for (int j = 0; j < n; ++j) { double s = std::pow(10.0, 150.0 - 300.0 * j / std::max(1, n - 1)); for (int i = 0; i < m; ++i) A[i + (size_t)j * m] *= s; }
  } else if (type == 5) {   // zero columns (tau = 0 paths) and real columns (alpha.i = 0)
    for (int j = 0; j < n; ++j) { if (j % 5 == 0) for (int i = 0; i < m; ++i) A[i + (size_t)j * m] = 0; else if (j % 5 == 1) for (int i = 0; i < m; ++i) A[i + (size_t)j * m] = cd(A[i + (size_t)j * m].real(), 0); }
  }
}
// device -> host through a plain memcpy; the V2 adapter completes it some ms after wait(): sleep (p4 STATUS)
static void d2h(sycl::queue &q, void *dst, const void *src, size_t bytes) { q.memcpy(dst, src, bytes).wait(); std::this_thread::sleep_for(std::chrono::milliseconds(30)); }
static void h2d(sycl::queue &q, void *dst, const void *src, size_t bytes) { q.memcpy(dst, src, bytes).wait(); }
static double maxabs(const std::vector<cd> &a) { double m = 0; for (auto &x : a) m = std::max(m, std::abs(x)); return m; }
static double reldiff(const std::vector<cd> &a, const std::vector<cd> &ref) { double d = 0; for (size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - ref[i])); double r = maxabs(ref); return r > 0 ? d / r : d; }

// true device time of f(0..reps-1) submitted back to back: closed by a kernel that sets a host-USM flag, the host spins on it
template <class F> static double chain_ms(sycl::queue &q, int reps, F &&f) {
  volatile int *flag = sycl::malloc_host<int>(1, q);
  *flag = 0;
  int *fl = (int *)flag;
  double t0 = now();
  for (int r = 0; r < reps; ++r) f(r);
  q.single_task([=] { *fl = 1; });
  while (*flag == 0) { }
  double t = now() - t0;
  sycl::free(fl, q);
  return 1e3 * t / reps;
}

static int nfail = 0;
static int plant = 0;

// ---- geqrf: vs host zgeqrf ----
static void test_geqrf(sycl::queue &q, int m, int n, int type, int nb) {
  const int k = std::min(m, n);
  std::vector<cd> A0; gen(type, m, n, A0, 12345 + m + n);
  std::vector<cd> H(A0), th(k);
  { MKL_INT info, mm = m, nn = n, lw = -1; cd wq; zgeqrf(&mm, &nn, H.data(), &mm, th.data(), &wq, &lw, &info); lw = (MKL_INT)wq.real(); std::vector<cd> w(lw); zgeqrf(&mm, &nn, H.data(), &mm, th.data(), w.data(), &lw, &info); }
  cd *dA = sycl::malloc_device<cd>((size_t)m * n, q), *dt = sycl::malloc_device<cd>(k, q);
  h2d(q, dA, A0.data(), sizeof(cd) * (size_t)m * n);
  int rc = alfdevqr::geqrf2(q, m, n, dA, m, dt, nb);
  q.wait_and_throw();
  std::vector<cd> D((size_t)m * n), td(k);
  d2h(q, D.data(), dA, sizeof(cd) * (size_t)m * n); d2h(q, td.data(), dt, sizeof(cd) * k);
  if (plant == 1) td[0] *= (1 + 1e-12);
  const bool skipLast = (m <= n) && std::abs(th[k - 1]) == 0 && std::abs(td[k - 1]) != 0;   // host MKL leaves a 1-row last reflector as identity
  double dR = 0, dV = 0, dT = 0;
  for (int j = 0; j < n; ++j) {
    double nr = 0, nv = 0, er = 0, ev = 0;
    for (int i = 0; i < m; ++i) {
      cd h = H[i + (size_t)j * m], d = D[i + (size_t)j * m];
      if (skipLast && i == k - 1) continue;
      if (i <= j) { nr = std::max(nr, std::abs(h)); er = std::max(er, std::abs(h - d)); }
      else { nv = std::max(nv, std::abs(h)); ev = std::max(ev, std::abs(h - d)); }
    }
    if (nr > 0) dR = std::max(dR, er / nr);
    if (nv > 0) dV = std::max(dV, ev / nv);
  }
  for (int j = 0; j < k - (skipLast ? 1 : 0); ++j) if (std::abs(th[j]) > 0) dT = std::max(dT, std::abs(th[j] - td[j]) / std::abs(th[j]));
  // orthogonality of the device Q: Q = ungqr(host) of the DEVICE reflectors
  auto orth = [&](const std::vector<cd> &F, const std::vector<cd> &tau) {
    MKL_INT kk = k, info, mm = m, ncol = m >= n ? n : m;
    std::vector<cd> Qm((size_t)m * ncol), t(tau);
    for (int j = 0; j < ncol; ++j) for (int i = 0; i < m; ++i) Qm[i + (size_t)j * m] = F[i + (size_t)j * m];
    cd wq; MKL_INT lw = -1; zungqr(&mm, &ncol, &kk, Qm.data(), &mm, t.data(), &wq, &lw, &info); lw = (MKL_INT)wq.real();
    std::vector<cd> w(lw); zungqr(&mm, &ncol, &kk, Qm.data(), &mm, t.data(), w.data(), &lw, &info);
    std::vector<cd> G((size_t)ncol * ncol); cd one(1, 0), zero(0, 0);
    cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, ncol, ncol, m, &one, Qm.data(), m, Qm.data(), m, &zero, G.data(), ncol);
    double mx = 0; for (int j = 0; j < ncol; ++j) for (int i = 0; i < ncol; ++i) mx = std::max(mx, std::abs(G[i + (size_t)j * ncol] - (i == j ? 1.0 : 0.0)));
    return mx;
  };
  const double oh = orth(H, th), od = orth(D, td);
  const bool ok = rc == 0 && dR <= 1e-13 && dV <= 1e-13 && dT <= 1e-13 && od <= 2 * oh + 1e-15;
  if (!ok) ++nfail;
  std::printf("%s geqrf m=%4d n=%4d type=%d nb=%2d rc=%d  dR=%.2e dV=%.2e dtau=%.2e  orth host=%.2e dev=%.2e%s\n", ok ? "PASS" : "FAIL", m, n, type, nb, rc, dR, dV, dT, oh, od, skipLast ? " [last row skipped]" : "");
  sycl::free(dA, q); sycl::free(dt, q);
}

// ---- unmqr: device apply of the HOST reflectors vs host zunmqr ----
static void test_unmqr(sycl::queue &q, char side, char trans, int m, int n, int k, int type) {
  const int nq = side == 'L' ? m : n;
  std::vector<cd> F; gen(type, nq, k, F, 777 + m + n + k);   // matrix to factor
  std::vector<cd> tau(k);
  { MKL_INT info, mm = nq, nn = k, lw = -1; cd wq; zgeqrf(&mm, &nn, F.data(), &mm, tau.data(), &wq, &lw, &info); lw = (MKL_INT)wq.real(); std::vector<cd> w(lw); zgeqrf(&mm, &nn, F.data(), &mm, tau.data(), w.data(), &lw, &info); }
  std::vector<cd> C0; gen(0, m, n, C0, 99 + m * 3 + n);
  std::vector<cd> R(C0);
  { MKL_INT info, mm = m, nn = n, kk = k, lda = nq, ldc = m, lw = -1; cd wq; zunmqr(&side, &trans, &mm, &nn, &kk, F.data(), &lda, tau.data(), R.data(), &ldc, &wq, &lw, &info); lw = (MKL_INT)wq.real(); std::vector<cd> w(lw);
    zunmqr(&side, &trans, &mm, &nn, &kk, F.data(), &lda, tau.data(), R.data(), &ldc, w.data(), &lw, &info); }
  std::vector<cd> tauin(tau);
  if (plant == 3) tauin[0] *= (1 + 1e-12);
  cd *dF = sycl::malloc_device<cd>(F.size(), q), *dt = sycl::malloc_device<cd>(k, q), *dC = sycl::malloc_device<cd>(C0.size(), q);
  h2d(q, dF, F.data(), sizeof(cd) * F.size()); h2d(q, dt, tauin.data(), sizeof(cd) * k); h2d(q, dC, C0.data(), sizeof(cd) * C0.size());
  int tr = trans == 'N' ? 0 : 2;
  if (plant == 2) tr = 2 - tr;
  int rc = alfdevqr::unmqr2(q, side == 'L' ? 0 : 1, tr, m, n, k, dF, nq, dt, dC, m);
  q.wait_and_throw();
  std::vector<cd> D(C0.size());
  d2h(q, D.data(), dC, sizeof(cd) * C0.size());
  const double d = reldiff(D, R);
  const bool ok = rc == 0 && d <= 1e-13;
  if (!ok) ++nfail;
  std::printf("%s unmqr %c%c m=%4d n=%4d k=%4d type=%d rc=%d  reldiff vs host=%.2e\n", ok ? "PASS" : "FAIL", side, trans, m, n, k, type, rc, d);
  sycl::free(dF, q); sycl::free(dt, q); sycl::free(dC, q);
}
// Q^H Q = I with the device apply only: apply Q (L,N) to the identity, then Q^H (L,C): result must be I; and Q from L,N is unitary
static void test_unitary(sycl::queue &q, int n) {
  std::vector<cd> A0; gen(0, n, n, A0, 5 + n);
  cd *dA = sycl::malloc_device<cd>((size_t)n * n, q), *dt = sycl::malloc_device<cd>(n, q), *dQ = sycl::malloc_device<cd>((size_t)n * n, q);
  h2d(q, dA, A0.data(), sizeof(cd) * (size_t)n * n);
  alfdevqr::geqrf2(q, n, n, dA, n, dt, 16);
  std::vector<cd> I((size_t)n * n, 0); for (int i = 0; i < n; ++i) I[i + (size_t)i * n] = 1;
  h2d(q, dQ, I.data(), sizeof(cd) * I.size());
  alfdevqr::unmqr2(q, 0, 0, n, n, n, dA, n, dt, dQ, n);   // Q
  q.wait_and_throw();
  std::vector<cd> Q(I.size()); d2h(q, Q.data(), dQ, sizeof(cd) * Q.size());
  std::vector<cd> G(I.size()); cd one(1, 0), zero(0, 0);
  cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, n, n, &one, Q.data(), n, Q.data(), n, &zero, G.data(), n);
  double mx = 0; for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) mx = std::max(mx, std::abs(G[i + (size_t)j * n] - (i == j ? 1.0 : 0.0)));
  // Q R = A
  std::vector<cd> D((size_t)n * n); d2h(q, D.data(), dA, sizeof(cd) * D.size());
  std::vector<cd> Rm((size_t)n * n, 0); for (int j = 0; j < n; ++j) for (int i = 0; i <= j; ++i) Rm[i + (size_t)j * n] = D[i + (size_t)j * n];
  std::vector<cd> QR(I.size()); cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Q.data(), n, Rm.data(), n, &zero, QR.data(), n);
  double dqr = reldiff(QR, A0);
  const bool ok = mx <= 1e-13 && dqr <= 1e-13;
  if (!ok) ++nfail;
  std::printf("%s unitary n=%4d  ||Q^H Q - I||=%.2e  ||QR-A||/||A||=%.2e\n", ok ? "PASS" : "FAIL", n, mx, dqr);
  sycl::free(dA, q); sycl::free(dt, q); sycl::free(dQ, q);
}

int main(int argc, char **argv) {
  std::string mode = argc > 1 ? argv[1] : "correct";
  int nb = argc > 2 ? atoi(argv[2]) : 16;
  plant = std::getenv("DEVQR_PLANT") ? atoi(std::getenv("DEVQR_PLANT")) : 0;
  sycl::queue q{sycl::gpu_selector_v, sycl::property::queue::in_order()};
  std::printf("device %s mode=%s nb=%d plant=%d\n", q.get_device().get_info<sycl::info::device::name>().c_str(), mode.c_str(), nb, plant);
  if (mode == "correct") {
    for (int type : {0, 1, 2, 4, 5}) for (int n : {36, 144, 576, 1024}) for (int b : {nb}) test_geqrf(q, n, n, type, b);
    for (int type : {0, 5}) { test_geqrf(q, 576, 144, type, nb); test_geqrf(q, 144, 576, type, nb); test_geqrf(q, 1000, 37, type, nb); test_geqrf(q, 37, 1000, type, nb); test_geqrf(q, 1300, 90, type, nb); test_geqrf(q, 1536, 64, type, nb); }
    for (int b : {8, 32}) { test_geqrf(q, 576, 576, 0, b); test_geqrf(q, 1024, 1024, 0, b); }
    test_geqrf(q, 5, 5, 0, nb); test_geqrf(q, 1, 1, 0, nb); test_geqrf(q, 17, 33, 2, nb); test_geqrf(q, 577, 577, 0, nb); test_geqrf(q, 130, 130, 1, nb);
    for (char side : {'L', 'R'}) for (char trans : {'N', 'C'}) {
      for (int n : {36, 144, 576, 1024}) test_unmqr(q, side, trans, n, n, n, 0);
      test_unmqr(q, side, trans, 576, 576, 576, 1); test_unmqr(q, side, trans, 576, 576, 576, 5);
      test_unmqr(q, side, trans, 200, 90, 70, 0); test_unmqr(q, side, trans, 90, 200, 60, 0); test_unmqr(q, side, trans, 1000, 33, 33, 0); test_unmqr(q, side, trans, 33, 1000, 33, 0);
      test_unmqr(q, side, trans, 1300, 100, 100, 0); test_unmqr(q, side, trans, 100, 1300, 100, 0); test_unmqr(q, side, trans, 7, 5, 5, 0); test_unmqr(q, side, trans, 1, 1, 1, 0); test_unmqr(q, side, trans, 577, 577, 300, 0);
    }
    for (int n : {64, 576, 1024}) test_unitary(q, n);
    std::printf("correctness: %d failed\n", nfail);
    return nfail ? 1 : 0;
  }
  if (mode == "time") {
    const int reps = 20;
    for (int n : {576, 1024}) {
      std::vector<cd> A0, C0; gen(0, n, n, A0, 7); gen(0, n, n, C0, 8);
      cd *A = sycl::malloc_device<cd>((size_t)n * n, q), *A00 = sycl::malloc_device<cd>((size_t)n * n, q), *t = sycl::malloc_device<cd>(n, q), *C = sycl::malloc_device<cd>((size_t)n * n, q), *C00 = sycl::malloc_device<cd>((size_t)n * n, q);
      int64_t ws = alfdevqr::scratch_size(n, n, 16), *dw0 = nullptr; (void)dw0;
      cd *w = sycl::malloc_device<cd>(ws, q);
      h2d(q, A00, A0.data(), sizeof(cd) * (size_t)n * n); h2d(q, C00, C0.data(), sizeof(cd) * (size_t)n * n);
      const size_t bytes = sizeof(cd) * (size_t)n * n;
      for (int rep = 0; rep < 2; ++rep) {   // first pass: JIT / warm-up
        double tc = chain_ms(q, reps, [&](int) { q.memcpy(A, A00, bytes); });
        double tg1 = chain_ms(q, reps, [&](int) { q.memcpy(A, A00, bytes); alfdevqr::geqrf(q, n, n, A, n, t, w, ws, 16, nullptr); }) - tc;
        double tg2 = chain_ms(q, reps, [&](int) { q.memcpy(A, A00, bytes); alfdevqr::geqrf2(q, n, n, A, n, t, nb); }) - tc;
        // factor once with geqrf2, then time unmqr on a copy of C
        q.memcpy(A, A00, bytes); alfdevqr::geqrf2(q, n, n, A, n, t, nb); q.wait();
        double tcc = chain_ms(q, reps, [&](int) { q.memcpy(C, C00, bytes); });
        double tu[4];
        int idx = 0;
        for (int side : {0, 1}) for (int trans : {0, 2}) tu[idx++] = chain_ms(q, reps, [&](int) { q.memcpy(C, C00, bytes); alfdevqr::unmqr2(q, side, trans, n, n, n, A, n, t, C, n); }) - tcc;
        if (rep) std::printf("TIME N=%4d  copy %.3f ms | geqrf old(devqr1) %.3f ms  new(geqrf2 nb=%d) %.3f ms | unmqr new L,N %.3f  L,C %.3f  R,N %.3f  R,C %.3f ms\n", n, tc, tg1, nb, tg2, tu[0], tu[1], tu[2], tu[3]);
      }
      sycl::free(A, q); sycl::free(A00, q); sycl::free(t, q); sycl::free(C, q); sycl::free(C00, q); sycl::free(w, q);
    }
    return 0;
  }
  if (mode == "apply") {   // apply kernel: time vs number of reflectors k, columns nc, rows M (fixed cost and cost per reflector)
    const int n = 1536;
    std::vector<cd> A0; gen(0, n, 128, A0, 3);
    cd *V = sycl::malloc_device<cd>((size_t)n * 128, q), *t = sycl::malloc_device<cd>(128, q), *Cc = sycl::malloc_device<cd>((size_t)n * 1024, q);
    h2d(q, V, A0.data(), sizeof(cd) * (size_t)n * 128);
    std::vector<cd> tt(128, cd(1.0, 0.1)); h2d(q, t, tt.data(), sizeof(cd) * 128);
    q.memset(Cc, 0, sizeof(cd) * (size_t)n * 1024).wait();
    for (int rep = 0; rep < 2; ++rep)
      for (int M : {128, 256, 384, 576, 1024}) for (int nc : {64, 512}) {
        double tk[3]; int ks[3] = {1, 16, 64};
        for (int i = 0; i < 3; ++i) tk[i] = chain_ms(q, 20, [&](int) { alfdevqr::apply_left<false>(q, M, ks[i], nc, true, true, (const double *)V, n, (const double *)t, (double *)Cc, n); });
        if (rep) std::printf("APPLY M=%4d nc=%3d: k=1 %7.1f us  k=16 %7.1f us  k=64 %7.1f us  -> %.2f us per reflector, fixed %.1f us\n", M, nc, 1e3 * tk[0], 1e3 * tk[1], 1e3 * tk[2], 1e3 * (tk[2] - tk[1]) / 48, 1e3 * (tk[1] - 16 * (tk[2] - tk[1]) / 48));
      }
    return 0;
  }
  if (mode == "panel") {
    const int n = 1024;
    std::vector<cd> A0; gen(0, n, 64, A0, 3);
    cd *P = sycl::malloc_device<cd>((size_t)n * 64, q), *P0 = sycl::malloc_device<cd>((size_t)n * 64, q), *t = sycl::malloc_device<cd>(64, q);
    h2d(q, P0, A0.data(), sizeof(cd) * (size_t)n * 64);
    for (int rep = 0; rep < 2; ++rep)
      for (int M : {64, 128, 256, 384, 512, 576, 768, 1024}) for (int jb : {8, 16, 32}) {
        if (jb > alfdevqr::panel2_maxjb(M)) continue;
        const size_t bytes = sizeof(cd) * (size_t)M * jb;
        double tc = chain_ms(q, 20, [&](int) { q.memcpy(P, P0, bytes); });
        double tp = chain_ms(q, 20, [&](int) { q.memcpy(P, P0, bytes); alfdevqr::panel2(q, M, jb, (double *)P, M, (double *)t); }) - tc;
        // apply of the jb reflectors to 512 columns
        cd *Cc = sycl::malloc_device<cd>((size_t)M * 512, q);
        q.memset(Cc, 0, sizeof(cd) * (size_t)M * 512).wait();
        double ta = chain_ms(q, 20, [&](int) { alfdevqr::apply_left<false>(q, M, jb, 512, true, true, (const double *)P, M, (const double *)t, (double *)Cc, M); });
        sycl::free(Cc, q);
        if (rep) std::printf("PANEL M=%4d jb=%2d: %8.1f us (%.2f us/col)   apply(512 cols) %8.1f us\n", M, jb, 1e3 * tp, 1e3 * tp / jb, 1e3 * ta);
      }
    return 0;
  }
  return 0;
}
