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
#include <sycl/ext/intel/experimental/grf_size_properties.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <cstdio>

namespace alfdevqr {
using cd = std::complex<double>;

struct Z { double r, i; };
inline Z zmul(Z a, Z b) { return {a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r}; }
inline Z zmulc(Z a, Z b) { return {a.r * b.r + a.i * b.i, a.r * b.i - a.i * b.r}; }   // conj(a) * b

constexpr int SGS = 16;                // sub-group size (lanes)
constexpr int MAX_NB = 48;             // jb sub-groups of 16 lanes = 768 work-items (64 -> 1024 work-items fails to launch on PVC)
constexpr int64_t MAX_M = 1536;        // panel height limit (SLM: 2 reflector buffers of M complex)

// dlapy3 of LAPACK
inline double dlapy3(double x, double y, double z) {
  const double hugeval = 1.7976931348623157e308;
  double xa = std::fabs(x), ya = std::fabs(y), za = std::fabs(z);
  double w = std::fmax(xa, std::fmax(ya, za));
  if (w == 0.0 || w > hugeval) return xa + ya + za;
  double a = xa / w, b = ya / w, c = za / w;
  return w * std::sqrt(a * a + b * b + c * c);
}
// sqrt(x^2+y^2+z^2) like dlapy3 (which scales by the largest |component|, three divisions): direct formula when the squares cannot
// overflow/underflow (differs from dlapy3 by rounding only), dlapy3 otherwise
inline double hypot3(double x, double y, double z) {
  double w = std::fmax(std::fabs(x), std::fmax(std::fabs(y), std::fabs(z)));
  if (w > 1e-140 && w < 1e140) return std::sqrt(x * x + y * y + z * z);
  return dlapy3(x, y, z);
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

// ---- panel kernel (register resident) -------------------------------------------------------------------------------
// Panel of M rows x jb columns at P (column-major, ld lda): ZGEQR2 + ZLARFT. Writes R/reflectors into P, tau(0:jb-1) and, if want_t,
// the explicit V (M x jb, ld M: zeros above the diagonal, ones on it) and the upper-triangular T (jb x jb, ld jb, zeros below).
// One work-group; every panel column lives in the REGISTERS of G sub-groups (16 lanes each, large GRF): lane l of sub-group (c,g)
// holds the rows (k*G+g)*16+l, k < KM, of column c for the whole factorisation (global memory latency, ~1 us, would otherwise
// dominate every column step; the kernel is bound by the instruction issue of the single Xe core it runs on, so the row loops
// are branch- and mask-free: the owner publishes the reflector zero-padded in SLM).
// Step i: the OWNER sub-group(s) of column i make the reflector (ZLARFG; the column and its norm were published by the sub-group(s)
// during step i-1: look-ahead), write the final column and the zero-padded reflector v (SLM) and tau; one work-group barrier; the
// sub-groups of the columns to its right apply H(i)^H (dot product by lane reduction [+ SLM exchange between the G sub-groups of a
// column], rank-1 update in registers); the next owner publishes its updated column's norm; the idle sub-groups l < i form
// Z(l,i) = v_l^H v_i from their registers, and sub-group 0 turns Z into the ZLARFT recurrence T(:,i) one step later.
template <int G, int KM>
inline void panel_kernel_t(sycl::queue &q, int M, int jb, double *P, int64_t lda, double *tau, double *Vx, double *Tm, bool want_t) {
  const int WI = jb * G * SGS;
  constexpr int MP = KM * G * SGS;   // padded column length
  q.submit([&](sycl::handler &h) {
    sycl::local_accessor<double, 1> vb(sycl::range<1>(2 * 2 * (size_t)MP), h);   // 2 buffers: zero-padded reflector v_i
    sycl::local_accessor<double, 1> tb(sycl::range<1>(4), h);                    // 2 buffers: tau_i
    sycl::local_accessor<double, 1> sx(sycl::range<1>(2 * G + 2), h);            // published: partial sums of squares [g], alpha (re,im)
    sycl::local_accessor<double, 1> wp(sycl::range<1>(2 * (size_t)jb * G + 8), h);   // partial dot products (G > 1)
    sycl::local_accessor<double, 1> zs(sycl::range<1>(2 * 2 * (size_t)G * jb), h);   // Z(l,i) partials [parity][g][l]
    sycl::local_accessor<double, 1> ts(sycl::range<1>(2 * (size_t)jb * jb), h);  // T
    sycl::ext::oneapi::experimental::properties props{sycl::ext::intel::experimental::grf_size<256>};
    h.parallel_for(sycl::nd_range<1>(WI, WI), props, [=](sycl::nd_item<1> it) [[sycl::reqd_sub_group_size(SGS)]] {
      auto sg = it.get_sub_group();
      const int sid = (int)sg.get_group_id()[0];
      const int c = sid / G, g = sid % G;
      const int lane = (int)sg.get_local_id()[0];
      const double safmin = 2.2250738585072014e-308 / 1.1102230246251565e-16;
      const double rsafmn = 1.0 / safmin;
      Z a[KM];
      {
        const double *src = P + 2 * (int64_t)c * lda;
#pragma unroll
        for (int k = 0; k < KM; ++k) { int r = (k * G + g) * SGS + lane; a[k] = r < M ? Z{src[2 * r], src[2 * r + 1]} : Z{0.0, 0.0}; }
      }
      Z my_scalx = {1.0, 0.0};
      // publish column c as the next owner: sx[2g], sx[2g+1] = (mx, s2) with partial sum of squares of rows > row = mx^2 * s2
      // (mx = 1 normally; the scaled form when the plain sum would leave [1e-280, 1e280]), and alpha = a[row]
      auto publish = [&](int row) {
        double s = 0.0;
#pragma unroll
        for (int k = 0; k < KM; ++k) {
          const int r = (k * G + g) * SGS + lane;
          s += (r > row && r < M) ? a[k].r * a[k].r + a[k].i * a[k].i : 0.0;
          if (r == row) { sx[2 * G] = a[k].r; sx[2 * G + 1] = a[k].i; }
        }
        s = sycl::reduce_over_group(sg, s, sycl::plus<double>());
        double mxg = 1.0;
        if (s != 0.0 && !(s > 1e-280 && s < 1e280)) {   // rare: scaled sum of squares of this sub-group's rows
          double mx = 0.0;
#pragma unroll
          for (int k = 0; k < KM; ++k) { const int r = (k * G + g) * SGS + lane; if (r > row && r < M) mx = std::fmax(mx, std::fmax(std::fabs(a[k].r), std::fabs(a[k].i))); }
          mx = sycl::reduce_over_group(sg, mx, sycl::maximum<double>());
          if (mx <= 1.7976931348623157e308) {
            double s2 = 0.0;
#pragma unroll
            for (int k = 0; k < KM; ++k) { const int r = (k * G + g) * SGS + lane; if (r > row && r < M) { double xr = a[k].r / mx, xi = a[k].i / mx; s2 += xr * xr + xi * xi; } }
            s = sycl::reduce_over_group(sg, s2, sycl::plus<double>());
            mxg = mx;
          } else { mxg = mx; s = 1.0; }   // inf / nan: propagate
        }
        if (lane == 0) { sx[2 * g] = mxg; sx[2 * g + 1] = s; }
      };
      for (int e = (int)it.get_local_linear_id(); e < jb * jb; e += WI) { ts[2 * e] = 0.0; ts[2 * e + 1] = 0.0; }
      if (c == 0) publish(0);
      sycl::group_barrier(it.get_group());

      // T column i from Z(:,i) (step i's partials) and tau_i: sub-group 0
      auto build_t = [&](int i) {
        const double ti = tb[2 * (i & 1)], tj = tb[2 * (i & 1) + 1];
        const double *z = &zs[2 * (i & 1) * G * jb];
        for (int l = lane; l <= i; l += SGS) {
          Z res = {0.0, 0.0};
          if (l == i) res = {ti, tj};
          else if (ti != 0.0 || tj != 0.0) {
            Z y = {0.0, 0.0};
            for (int p = l; p < i; ++p) {
              Z zz = {0.0, 0.0};
              for (int gg = 0; gg < G; ++gg) { zz.r += z[2 * (gg * jb + p)]; zz.i += z[2 * (gg * jb + p) + 1]; }
              Z pp = zmul({ts[2 * (l + p * jb)], ts[2 * (l + p * jb) + 1]}, zz);
              y.r += pp.r; y.i += pp.i;
            }
            res = zmul({-ti, -tj}, y);
          }
          ts[2 * (l + i * jb)] = res.r; ts[2 * (l + i * jb) + 1] = res.i;
        }
        sycl::group_barrier(sg);
      };

      for (int i = 0; i < jb; ++i) {
        const int buf = i & 1;
        double *v = &vb[buf * 2 * MP];
        if (G > 1) sycl::group_barrier(it.get_group());   // published sums/alpha visible to both sub-groups of the owner
        if (c == i) {   // ---- owner: ZLARFG, final column, zero-padded reflector
          if (G == 1) sycl::group_barrier(sg);
          Z alpha = {sx[2 * G], sx[2 * G + 1]};
          double xn;
          {
            double mxs = 0.0;
            bool scaled = false;
#pragma unroll
            for (int gg = 0; gg < G; ++gg) { mxs = std::fmax(mxs, sx[2 * gg]); scaled = scaled || sx[2 * gg] != 1.0; }
            if (!scaled) {
              double s_tot = 0.0;
#pragma unroll
              for (int gg = 0; gg < G; ++gg) s_tot += sx[2 * gg + 1];
              xn = s_tot == 0.0 ? 0.0 : std::sqrt(s_tot);
            } else if (!(mxs <= 1.7976931348623157e308)) xn = mxs;
            else {
              double ss = 0.0;
#pragma unroll
              for (int gg = 0; gg < G; ++gg) { double q2 = sx[2 * gg] / mxs; ss += q2 * q2 * sx[2 * gg + 1]; }
              xn = mxs * std::sqrt(ss);
            }
          }
          Z t = {0.0, 0.0}, scalx = {1.0, 0.0};
          double beta = alpha.r;
          const bool gen = !(xn == 0.0 && alpha.i == 0.0);
          if (gen) {
            beta = -std::copysign(hypot3(alpha.r, alpha.i, xn), alpha.r);
            int knt = 0;
            double skn = 1.0;
            if (std::fabs(beta) < safmin) {   // LAPACK: scale x, alpha up until beta is safely representable, then recompute
              do { ++knt; skn *= rsafmn; beta *= rsafmn; } while (std::fabs(beta) < safmin && knt < 20);
              alpha.r *= skn; alpha.i *= skn; xn *= skn;
              beta = -std::copysign(hypot3(alpha.r, alpha.i, xn), alpha.r);
            }
            t = {(beta - alpha.r) / beta, -alpha.i / beta};
            Z sc1 = zinv({alpha.r - beta, alpha.i});
            scalx = {sc1.r * skn, sc1.i * skn};
            for (int k2 = 0; k2 < knt; ++k2) beta *= safmin;
          }
          my_scalx = scalx;
          double *dst = P + 2 * (int64_t)c * lda;
          double *vxc = Vx + 2 * (int64_t)c * M;
#pragma unroll
          for (int k = 0; k < KM; ++k) {
            const int r = (k * G + g) * SGS + lane;
            const Z vv = zmul(a[k], scalx);
            const Z o = r < i ? a[k] : vv;
            if (r < M && r != i) { dst[2 * r] = o.r; dst[2 * r + 1] = o.i; }
            const Z vz = r < i ? Z{0.0, 0.0} : (r == i ? Z{1.0, 0.0} : vv);
            v[2 * r] = vz.r; v[2 * r + 1] = vz.i;
            if (want_t && r < M) { vxc[2 * r] = vz.r; vxc[2 * r + 1] = vz.i; }
          }
          if (g == 0 && lane == 0) { dst[2 * i] = beta; dst[2 * i + 1] = 0.0; tau[2 * i] = t.r; tau[2 * i + 1] = t.i; tb[2 * buf] = t.r; tb[2 * buf + 1] = t.i; }
        }
        sycl::group_barrier(it.get_group());
        Z w = {0.0, 0.0};
        double tr = 0.0, ti = 0.0;
        bool gen = false;
        if (c > i) {
          tr = tb[2 * buf]; ti = tb[2 * buf + 1];
          gen = (tr != 0.0 || ti != 0.0);
          if (gen) {   // pass A: w = v^H a
#pragma unroll
            for (int k = 0; k < KM; ++k) {
              const int r = (k * G + g) * SGS + lane;
              const Z vk = {v[2 * r], v[2 * r + 1]};
              const Z p = zmulc(vk, a[k]);
              w.r += p.r; w.i += p.i;
            }
            w.r = sycl::reduce_over_group(sg, w.r, sycl::plus<double>());
            w.i = sycl::reduce_over_group(sg, w.i, sycl::plus<double>());
          }
        }
        if (G > 1) {   // the G sub-groups of a column combine their partial dot products
          if (gen && lane == 0) { wp[2 * sid] = w.r; wp[2 * sid + 1] = w.i; }
          sycl::group_barrier(it.get_group());
          if (gen) {
            w = {0.0, 0.0};
#pragma unroll
            for (int gg = 0; gg < G; ++gg) { w.r += wp[2 * (c * G + gg)]; w.i += wp[2 * (c * G + gg) + 1]; }   // same order in all G sub-groups
          }
        }
        if (c > i) {
          if (gen) {   // pass B: a -= conj(tau) w v
            const Z f = zmul({tr, -ti}, w);
#pragma unroll
            for (int k = 0; k < KM; ++k) {
              const int r = (k * G + g) * SGS + lane;
              const Z vk = {v[2 * r], v[2 * r + 1]};
              const Z p = zmul(f, vk);
              a[k].r -= p.r; a[k].i -= p.i;
            }
          }
          if (c == i + 1) publish(i + 1);   // next owner
        } else if (c < i && want_t) {   // idle sub-groups: Z(c,i) = v_c^H v_i from the registers (v_c = a * scalx, rows >= i > c)
          Z z = {0.0, 0.0};
#pragma unroll
          for (int k = 0; k < KM; ++k) {
            const int r = (k * G + g) * SGS + lane;
            const Z vk = {v[2 * r], v[2 * r + 1]};
            const Z p = zmulc(a[k], vk);
            z.r += p.r; z.i += p.i;
          }
          z.r = sycl::reduce_over_group(sg, z.r, sycl::plus<double>());
          z.i = sycl::reduce_over_group(sg, z.i, sycl::plus<double>());
          z = zmulc(my_scalx, z);   // conj(scalx_c) * sum conj(a) v_i
          if (lane == 0) { zs[2 * ((i & 1) * G * jb + g * jb + c)] = z.r; zs[2 * ((i & 1) * G * jb + g * jb + c) + 1] = z.i; }
        }
        // sub-group 0 builds T(:, i-1) from the previous step's Z (visible since the barrier above)
        if (want_t && c == 0 && g == 0 && i >= 1) build_t(i - 1);
      }
      sycl::group_barrier(it.get_group());
      if (want_t && c == 0 && g == 0) {
        build_t(jb - 1);
        for (int e = lane; e < jb * jb; e += SGS) { Tm[2 * e] = ts[2 * e]; Tm[2 * e + 1] = ts[2 * e + 1]; }
      }
    });
  });
}

// configuration of a panel: widest panel for the height M (sub-groups per column G, register slots KM); 0 = not handled
inline int panel_maxjb(int64_t M) { return M <= 768 ? 32 : (M <= 1536 ? 16 : 0); }
inline void panel_kernel(sycl::queue &q, int M, int jb, double *P, int64_t lda, double *tau, double *Vx, double *Tm, bool want_t) {
  if (M <= 768) {
    if (M <= 256) panel_kernel_t<1, 16>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else if (M <= 384) panel_kernel_t<1, 24>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else if (M <= 512) panel_kernel_t<1, 32>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else if (M <= 640) panel_kernel_t<1, 40>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else panel_kernel_t<1, 48>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
  } else {
    if (M <= 1024) panel_kernel_t<2, 32>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else if (M <= 1280) panel_kernel_t<2, 40>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
    else panel_kernel_t<2, 48>(q, M, jb, P, lda, tau, Vx, Tm, want_t);
  }
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
  if (m > MAX_M || panel_maxjb(m) == 0 || lwork < scratch_size(m, n, nb)) return -1;
  static const bool prof = std::getenv("ALFDEVQR_PROF") != nullptr;   // per-stage host timing with a wait after every launch (diagnostic)
  static double tp[4] = {0, 0, 0, 0}; static long ncall = 0;
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  double tl = prof ? now() : 0.0;
  auto lap = [&](int k) { if (prof) { q.wait(); double t = now(); tp[k] += t - tl; tl = t; } };
  cd *Vx = work, *Tm = Vx + m * nb, *W = Tm + (int64_t)nb * nb, *W2 = W + (int64_t)nb * n;
  sycl::event ev;
  for (int64_t j = 0, jbp; j < k; j += jbp) {
    const int jb = (int)std::min<int64_t>(std::min(nb, panel_maxjb(m - j)), k - j);
    jbp = jb;
    const int M = (int)(m - j);
    const int64_t nc = n - j - jb;
    cd *P = A + j + j * lda;
    panel_kernel(q, M, jb, (double *)P, lda, (double *)(tau + j), (double *)Vx, (double *)Tm, nc > 0);
    lap(0);
    if (nc > 0) {
      cd *C = A + j + (j + jb) * lda;
      oneapi::mkl::blas::column_major::gemm(q, transpose::conjtrans, transpose::nontrans, jb, nc, M, cd(1.0), Vx, M, C, lda, cd(0.0), W, jb);
      lap(1);
      oneapi::mkl::blas::column_major::gemm(q, transpose::conjtrans, transpose::nontrans, jb, nc, jb, cd(1.0), Tm, jb, W, jb, cd(0.0), W2, jb);
      lap(2);
      ev = oneapi::mkl::blas::column_major::gemm(q, transpose::nontrans, transpose::nontrans, M, nc, jb, cd(-1.0), Vx, M, W2, jb, cd(1.0), C, lda);
      lap(3);
    }
  }
  if (prof) { if (++ncall % 50 == 0) std::fprintf(stderr, "ALFDEVQR_PROF calls %ld  s total: panel(+Z,T) %.3f gemm1 %.3f gemm2 %.3f gemm3 %.3f  (m=%lld n=%lld nb=%d)\n", ncall, tp[0], tp[1], tp[2], tp[3], (long long)m, (long long)n, nb); }
  if (last) *last = ev;
  return 0;
}
}  // namespace alfdevqr
