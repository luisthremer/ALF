// Device-only blocked Householder QR (complex double), the SAME algorithm as reference LAPACK ZGEQRF:
//   panel of nb columns : ZGEQR2 (ZLARFG + ZLARF, right-looking), ZLARFT (forward, columnwise)
//   trailing update     : ZLARFB('Left','Conjugate transpose','Forward','Columnwise') = C -= V (T^H (V^H C))
// Everything runs on the device queue: no host panel (oneMKL's ZGEQRF is a hybrid whose panels run on the host CPU).
// Output layout = LAPACK (R in the upper triangle, reflectors below the diagonal with unit diagonal implied, tau), so that
// ZUNMQR/ZUNGQR (oneMKL) and the k_* kernels of ALF consume it unchanged.
//
// Panel kernel: ONE work-group, one sub-group (16 lanes) per panel column; the owner of column i makes the reflector,
// publishes it in SLM, one work-group barrier per column, the sub-groups of the columns to its right apply it (dot product
// by sub-group reduction, rank-1 update). ZLARFT's V^H V is formed after the loop, T from the recurrence by one sub-group.
// The trailing update uses oneMKL device GEMMs on an explicit V (zeros above, ones on the diagonal).
//
// Used by alf_la_sycl_shim.cpp (ALFSH_DEVQR=1) and by the standalone test docs/supermuc_large_prompt/phase2/team/p4/devqr_test.cpp.
#pragma once
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cstdlib>

namespace alfdevqr {
using cd = std::complex<double>;

struct Z { double r, i; };
inline Z zmul(Z a, Z b) { return {a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r}; }
inline Z zmulc(Z a, Z b) { return {a.r * b.r + a.i * b.i, a.r * b.i - a.i * b.r}; }   // conj(a) * b

constexpr int SGS = 16;                // sub-group size (lanes)
constexpr int MAX_NB = 48;             // jb sub-groups of 16 lanes = 768 work-items (64 -> 1024 work-items fails to launch on PVC)
constexpr int64_t MAX_M = 4096;        // panel height limit (SLM: 2 reflector buffers of M complex)

// dlapy3 of LAPACK
inline double dlapy3(double x, double y, double z) {
  const double hugeval = 1.7976931348623157e308;
  double xa = std::fabs(x), ya = std::fabs(y), za = std::fabs(z);
  double w = std::fmax(xa, std::fmax(ya, za));
  if (w == 0.0 || w > hugeval) return xa + ya + za;
  double a = xa / w, b = ya / w, c = za / w;
  return w * std::sqrt(a * a + b * b + c * c);
}
// zladiv(1, y) (Smith)
inline Z zinv(Z y) {
  if (std::fabs(y.i) < std::fabs(y.r)) {
    double e = y.i / y.r, f = y.r + y.i * e;
    return {1.0 / f, -e / f};
  } else {
    double e = y.r / y.i, f = y.r * e + y.i;
    return {e / f, -1.0 / f};
  }
}

// Panel of M rows x jb columns at P (column-major, ld lda). Writes reflectors/R into P, tau(0:jb-1); if want_t also the explicit
// V (M x jb, ld M) and the upper-triangular T (jb x jb, ld jb, zeros below the diagonal).
inline void panel_kernel(sycl::queue &q, int M, int jb, double *P, int64_t lda, double *tau, double *Vx, double *Tm, bool want_t) {
  const int WI = jb * SGS;
  q.submit([&](sycl::handler &h) {
    sycl::local_accessor<double, 1> vb(sycl::range<1>(2 * 2 * (size_t)M), h);        // 2 reflector buffers (re,im)
    sycl::local_accessor<double, 1> scb(sycl::range<1>(8), h);                       // tau (re,im) per buffer
    sycl::local_accessor<double, 1> zt(sycl::range<1>(2 * (size_t)jb * jb), h);      // Z = V^H V (upper), later unused
    sycl::local_accessor<double, 1> tt(sycl::range<1>(2 * (size_t)jb * jb), h);      // T
    h.parallel_for(sycl::nd_range<1>(WI, WI), [=](sycl::nd_item<1> it) [[sycl::reqd_sub_group_size(SGS)]] {
      auto sg = it.get_sub_group();
      const int c = (int)sg.get_group_id()[0];
      const int lane = (int)sg.get_local_id()[0];
      const int kk = jb < M ? jb : M;
      double *a = P + 2 * (int64_t)c * lda;
      const double safmin = 2.2250738585072014e-308 / 1.1102230246251565e-16;
      const double rsafmn = 1.0 / safmin;
      for (int i = 0; i < kk; ++i) {
        const int buf = i & 1;
        double *v = &vb[buf * 2 * M];
        if (c == i) {
          // ---- ZLARFG on a(i:M)
          Z alpha = {a[2 * i], a[2 * i + 1]};
          double s = 0.0;
          for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
            if (r > i) { double xr = a[2 * r], xi = a[2 * r + 1]; s += xr * xr + xi * xi; }
          s = sycl::reduce_over_group(sg, s, sycl::plus<double>());
          double xnorm;
          if (s > 1e-280 && s < 1e280) xnorm = std::sqrt(s);
          else {   // scaled sum of squares (rare)
            double mx = 0.0;
            for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
              if (r > i) mx = std::fmax(mx, std::fmax(std::fabs(a[2 * r]), std::fabs(a[2 * r + 1])));
            mx = sycl::reduce_over_group(sg, mx, sycl::maximum<double>());
            if (mx == 0.0 || !(mx <= 1.7976931348623157e308)) xnorm = mx;
            else {
              double s2 = 0.0;
              for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
                if (r > i) { double xr = a[2 * r] / mx, xi = a[2 * r + 1] / mx; s2 += xr * xr + xi * xi; }
              s2 = sycl::reduce_over_group(sg, s2, sycl::plus<double>());
              xnorm = mx * std::sqrt(s2);
            }
          }
          Z t = {0.0, 0.0};
          Z scal = {1.0, 0.0};
          double beta = alpha.r;
          int knt = 0;
          bool gen = !(xnorm == 0.0 && alpha.i == 0.0);
          if (gen) {
            beta = -std::copysign(dlapy3(alpha.r, alpha.i, xnorm), alpha.r);
            if (std::fabs(beta) < safmin) {   // xnorm, beta may be inaccurate: scale x and recompute (LAPACK loop)
              do {
                ++knt;
                for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
                  if (r > i) { a[2 * r] *= rsafmn; a[2 * r + 1] *= rsafmn; }
                beta *= rsafmn; alpha.r *= rsafmn; alpha.i *= rsafmn;
              } while (std::fabs(beta) < safmin && knt < 20);
              double s3 = 0.0, mx = 0.0;
              for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
                if (r > i) { mx = std::fmax(mx, std::fmax(std::fabs(a[2 * r]), std::fabs(a[2 * r + 1]))); }
              mx = sycl::reduce_over_group(sg, mx, sycl::maximum<double>());
              if (mx > 0.0) {
                for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
                  if (r > i) { double xr = a[2 * r] / mx, xi = a[2 * r + 1] / mx; s3 += xr * xr + xi * xi; }
                s3 = sycl::reduce_over_group(sg, s3, sycl::plus<double>());
                xnorm = mx * std::sqrt(s3);
              } else xnorm = 0.0;
              beta = -std::copysign(dlapy3(alpha.r, alpha.i, xnorm), alpha.r);
            }
            t = {(beta - alpha.r) / beta, -alpha.i / beta};
            Z d = {alpha.r - beta, alpha.i};
            scal = zinv(d);
          }
          if (gen) {
            for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
              if (r > i) {
                Z x = zmul({a[2 * r], a[2 * r + 1]}, scal);
                a[2 * r] = x.r; a[2 * r + 1] = x.i;
                v[2 * r] = x.r; v[2 * r + 1] = x.i;
              }
            for (int k2 = 0; k2 < knt; ++k2) beta *= safmin;
          } else {   // H = I: x is exactly zero
            for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
              if (r > i) { v[2 * r] = a[2 * r]; v[2 * r + 1] = a[2 * r + 1]; }
          }
          if (lane == 0) {
            a[2 * i] = beta; a[2 * i + 1] = 0.0;
            v[2 * i] = 1.0; v[2 * i + 1] = 0.0;
            tau[2 * i] = t.r; tau[2 * i + 1] = t.i;
            scb[buf * 2] = t.r; scb[buf * 2 + 1] = t.i;
          }
        }
        sycl::group_barrier(it.get_group());
        if (c > i) {
          const double tr = scb[buf * 2], ti = scb[buf * 2 + 1];
          if (tr != 0.0 || ti != 0.0) {
            double wr = 0.0, wi = 0.0;
            for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
              if (r >= i) {
                Z vc = {v[2 * r], v[2 * r + 1]}, x = {a[2 * r], a[2 * r + 1]};
                Z p = zmulc(vc, x);
                wr += p.r; wi += p.i;
              }
            wr = sycl::reduce_over_group(sg, wr, sycl::plus<double>());
            wi = sycl::reduce_over_group(sg, wi, sycl::plus<double>());
            Z f = zmul({tr, -ti}, {wr, wi});   // conj(tau) * w
            for (int r = (i / SGS) * SGS + lane; r < M; r += SGS)
              if (r >= i) {
                Z p = zmul(f, {v[2 * r], v[2 * r + 1]});
                a[2 * r] -= p.r; a[2 * r + 1] -= p.i;
              }
          }
        }
      }
      if (!want_t) return;
      // explicit V (zeros above the diagonal, 1 on it) of column c
      double *vx = Vx + 2 * (int64_t)c * M;
      for (int r = lane; r < M; r += SGS) {
        double xr = 0.0, xi = 0.0;
        if (r == c) xr = 1.0;
        else if (r > c) { xr = a[2 * r]; xi = a[2 * r + 1]; }
        vx[2 * r] = xr; vx[2 * r + 1] = xi;
      }
      sycl::group_barrier(it.get_group());
      // Z(l,c) = V(:,l)^H V(:,c), l < c (rows >= c)
      for (int l = 0; l < c; ++l) {
        const double *al = Vx + 2 * (int64_t)l * M;
        double zr = 0.0, zi = 0.0;
        for (int r = (c / SGS) * SGS + lane; r < M; r += SGS)
          if (r >= c) {
            Z p = zmulc({al[2 * r], al[2 * r + 1]}, {vx[2 * r], vx[2 * r + 1]});
            zr += p.r; zi += p.i;
          }
        zr = sycl::reduce_over_group(sg, zr, sycl::plus<double>());
        zi = sycl::reduce_over_group(sg, zi, sycl::plus<double>());
        if (lane == 0) { zt[2 * (l + c * jb)] = zr; zt[2 * (l + c * jb) + 1] = zi; }
      }
      sycl::group_barrier(it.get_group());
      if (c == 0) {   // T recurrence (ZLARFT forward columnwise): T(0:i-1,i) = -tau_i T(0:i-1,0:i-1) Z(0:i-1,i)
        for (int i = 0; i < jb; ++i) {
          Z t = {tau[2 * i], tau[2 * i + 1]};
          for (int l = lane; l < jb; l += SGS) {
            Z res = {0.0, 0.0};
            if (l < i && (t.r != 0.0 || t.i != 0.0)) {
              Z y = {0.0, 0.0};
              for (int p = l; p < i; ++p) {
                Z pp = zmul({tt[2 * (l + p * jb)], tt[2 * (l + p * jb) + 1]}, {zt[2 * (p + i * jb)], zt[2 * (p + i * jb) + 1]});
                y.r += pp.r; y.i += pp.i;
              }
              Z m = zmul({-t.r, -t.i}, y);
              res = m;
            } else if (l == i) res = t;
            tt[2 * (l + i * jb)] = res.r; tt[2 * (l + i * jb) + 1] = res.i;
          }
          sycl::group_barrier(sg);
        }
        for (int e = lane; e < jb * jb; e += SGS) {
          int l = e % jb, cc = e / jb;
          Tm[2 * e] = l <= cc ? tt[2 * e] : 0.0;
          Tm[2 * e + 1] = l <= cc ? tt[2 * e + 1] : 0.0;
        }
      }
    });
  });
}

// naive reference kernels (DEVQR_NAIVE=1, debugging): C(m x n) = alpha*op(A) op(B) + beta*C, one work-item per output
inline void naive_gemm(sycl::queue &q, bool ca, int64_t m, int64_t n, int64_t k, double alr, const double *A, int64_t lda, const double *B, int64_t ldb, double bt, double *C, int64_t ldc) {
  q.parallel_for(sycl::range<2>(n, m), [=](sycl::id<2> id) {
    int64_t j = id[0], i = id[1];
    double sr = 0, si = 0;
    for (int64_t l = 0; l < k; ++l) {
      Z a = ca ? Z{A[2 * (l + i * lda)], -A[2 * (l + i * lda) + 1]} : Z{A[2 * (i + l * lda)], A[2 * (i + l * lda) + 1]};
      Z b = {B[2 * (l + j * ldb)], B[2 * (l + j * ldb) + 1]};
      Z p = zmul(a, b);
      sr += p.r; si += p.i;
    }
    C[2 * (i + j * ldc)] = bt * C[2 * (i + j * ldc)] + alr * sr;
    C[2 * (i + j * ldc) + 1] = bt * C[2 * (i + j * ldc) + 1] + alr * si;
  });
}

// scratch (complex<double> elements) needed by geqrf() for an m x n matrix with block size nb
inline int64_t scratch_size(int64_t m, int64_t n, int nb) { return m * nb + (int64_t)nb * nb + 2 * (int64_t)nb * n; }

// ZGEQRF on the device. Returns 0, or -1 (not handled: caller falls back: m too large for the panel kernel, scratch too small)
inline int geqrf(sycl::queue &q, int64_t m, int64_t n, cd *A, int64_t lda, cd *tau, cd *work, int64_t lwork, int nb,
                 sycl::event *last) {
  using oneapi::mkl::transpose;
  const int64_t k = std::min(m, n);
  if (k == 0) return 0;
  if (nb < 1) nb = 1;
  if (nb > MAX_NB) nb = MAX_NB;
  if (m > MAX_M || lwork < scratch_size(m, n, nb)) return -1;
  cd *Vx = work, *Tm = Vx + m * nb, *W = Tm + (int64_t)nb * nb, *W2 = W + (int64_t)nb * n;
  sycl::event ev;
  for (int64_t j = 0; j < k; j += nb) {
    const int jb = (int)std::min<int64_t>(nb, k - j);
    const int M = (int)(m - j);
    const int64_t nc = n - j - jb;
    cd *P = A + j + j * lda;
    static const bool dbg_sync = std::getenv("DEVQR_SYNC") != nullptr;
    panel_kernel(q, M, jb, (double *)P, lda, (double *)(tau + j), (double *)Vx, (double *)Tm, nc > 0);
    if (nc > 0) {
      cd *C = A + j + (j + jb) * lda;
      static const bool naive = std::getenv("DEVQR_NAIVE") != nullptr;
      if (naive) {
        naive_gemm(q, true, jb, nc, M, 1.0, (double *)Vx, M, (double *)C, lda, 0.0, (double *)W, jb);
        naive_gemm(q, true, jb, nc, jb, 1.0, (double *)Tm, jb, (double *)W, jb, 0.0, (double *)W2, jb);
        naive_gemm(q, false, M, nc, jb, -1.0, (double *)Vx, M, (double *)W2, jb, 1.0, (double *)C, lda);
        continue;
      }
      if (dbg_sync) q.wait();
      oneapi::mkl::blas::column_major::gemm(q, transpose::conjtrans, transpose::nontrans, jb, nc, M, cd(1.0), Vx, M, C, lda, cd(0.0), W, jb);
      if (dbg_sync) q.wait();
      oneapi::mkl::blas::column_major::gemm(q, transpose::conjtrans, transpose::nontrans, jb, nc, jb, cd(1.0), Tm, jb, W, jb, cd(0.0), W2, jb);
      if (dbg_sync) q.wait();
      ev = oneapi::mkl::blas::column_major::gemm(q, transpose::nontrans, transpose::nontrans, M, nc, jb, cd(-1.0), Vx, M, W2, jb, cd(1.0), C, lda);
    }
  }
  if (last) *last = ev;
  return 0;
}
}  // namespace alfdevqr
