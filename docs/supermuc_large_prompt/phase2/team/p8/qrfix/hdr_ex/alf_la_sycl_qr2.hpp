// Phase 5: T-free device QR and device ZUNMQR (complex double) for the Max 1550, after the HIP version of the AMD branch
// (alf_hip_qr.hip.cpp). Same LAPACK conventions as alf_la_sycl_devqr.hpp (ZLARFG semantics, R above the diagonal, reflectors below with unit
// diagonal implied, tau), so ZUNMQR/ZUNGQR of oneMKL and the callers in ALF consume the result unchanged.
//   panel kernel  (panel2_t) : ZGEQR2 only. ONE work-group, one sub-group (16 lanes) [or G sub-groups] per panel column, the column in registers,
//                              one work-group barrier per column. No ZLARFT: the phase-4 kernel spent most of its 7 us per column step in the
//                              T recurrence (sub-group 0, O(i) dependent SLM loads) and in the Z(c,i) dot products of the idle sub-groups.
//   apply kernel  (apply2_t) : applies a range of reflectors to independent vectors (columns of C, or conj-transposed rows of C for side R).
//                              One sub-group (Q row slots per lane) [G sub-groups, SLM exchange] per NV columns, the columns in registers,
//                              reflectors one after the other (dot by lane reduction, rank-1 update). No T factor, no GEMM, one launch for the
//                              trailing update of a QR panel and one for a whole ZUNMQR.
#pragma once
#include "alf_la_sycl_devqr.hpp"

namespace alfdevqr {
#ifdef QR2_TESTDELAY
inline int qr2_env_int(const char *n) { const char *e = std::getenv(n); return e ? std::atoi(e) : 0; }
#endif

// two-value lane reduction (xor butterfly: every lane gets the bitwise identical total)
inline void red2(sycl::sub_group &sg, double &a, double &b) {
#pragma unroll
  for (int m = 8; m >= 1; m >>= 1) {
    a += sycl::permute_group_by_xor(sg, a, m);
    b += sycl::permute_group_by_xor(sg, b, m);
  }
}

// ---- panel: ZGEQR2 of M x jb at P ----------------------------------------------------------------------------------
// lane l of sub-group (c,g) holds the rows (k*G+g)*16+l, k < KM, of column c. Step i: the owner (sub-group(s) of column i) makes the
// reflector from the norm published during step i-1 (look-ahead), writes it zero-padded to SLM, one work-group barrier, the sub-groups of the
// columns to its right apply H(i)^H (dot product by lane reduction [+ SLM exchange between the G sub-groups of a column], rank-1 update), the
// owner of step i+1 publishes the norm of its updated column. A finished column is written to global memory after the barrier (off the critical path).
template <int G, int KM>
inline void panel2_t(sycl::queue &q, int M, int jb, double *P, int64_t lda, double *tau) {
  const int WI = jb * G * SGS;
  constexpr int MP = KM * G * SGS;
  constexpr bool KEEPV = KM <= 20;
#ifdef QR2_TESTDELAY
  const int pdly_ = qr2_env_int("QR2_PDELAY");
#endif
  q.submit([&](sycl::handler &h) {
    sycl::local_accessor<sycl::double2, 1> vb(sycl::range<1>(2 * (size_t)MP), h);     // 2 buffers: zero-padded reflector v_i (complex, 16-byte loads)
    sycl::local_accessor<double, 1> tb(sycl::range<1>(4), h);                        // 2 buffers (i&1): tau_i
    sycl::local_accessor<double, 1> sx(sycl::range<1>(2 * G + 2), h);                // published: (mx, s2) per sub-group, alpha
    sycl::local_accessor<double, 1> wp(sycl::range<1>(2 * (size_t)jb * G + 8), h);   // partial dot products (G > 1)
    sycl::ext::oneapi::experimental::properties props{sycl::ext::intel::experimental::grf_size<256>};
    h.parallel_for(sycl::nd_range<1>(WI, WI), props, [=](sycl::nd_item<1> it) [[sycl::reqd_sub_group_size(SGS)]] {
      auto sg = it.get_sub_group();
      const int sid = (int)sg.get_group_id()[0];
      const int c = sid / G, g = sid % G;
      const int lane = (int)sg.get_local_id()[0];
      const double safmin = 2.2250738585072014e-308 / 1.1102230246251565e-16;
      const double rsafmn = 1.0 / safmin;
      Z a[KM];
      Z vreg[KEEPV ? KM : 1];
      {
        const double *src = P + 2 * (int64_t)c * lda;
#pragma unroll
        for (int k = 0; k < KM; ++k) { int r = (k * G + g) * SGS + lane; const int rc = r < M ? r : 0; a[k] = Z{src[2 * rc], src[2 * rc + 1]}; if (r >= M) a[k] = Z{0.0, 0.0}; }
      }
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
#ifndef ABL_B
        s = sycl::reduce_over_group(sg, s, sycl::plus<double>());
#endif
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
      if (c == 0) publish(0);
      it.barrier(sycl::access::fence_space::local_space);

      for (int i = 0; i < jb; ++i) {
        const int buf = i & 1;
        const int vo = buf * MP;   // direct accessor indexing: a pointer into SLM becomes a generic pointer (slow generic loads)
        if (G > 1) it.barrier(sycl::access::fence_space::local_space);   // published sums/alpha visible to all sub-groups of the owner
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
#ifdef ABL_A
          if (gen) { beta = alpha.r + xn; t = {0.5 * beta, 0.1}; scalx = {1.1 + beta, 0.2}; }
          if (false)
#else
          if (gen)
#endif
          {
            beta = -std::copysign(hypot3(alpha.r, alpha.i, xn), alpha.r);
            int knt = 0;
            double skn = 1.0;
            if (std::fabs(beta) < safmin) {   // LAPACK: scale x, alpha up until beta is safely representable, then recompute
              do { ++knt; skn *= rsafmn; beta *= rsafmn; } while (std::fabs(beta) < safmin && knt < 20);
              alpha.r *= skn; alpha.i *= skn; xn *= skn;
              beta = -std::copysign(hypot3(alpha.r, alpha.i, xn), alpha.r);
            }
            if (knt == 0 && std::fabs(beta) > 1e-100 && std::fabs(beta) < 1e100) {
              // fast path (no overflow/underflow possible: |alpha - beta| >= |beta|): the two divisions are independent (the LAPACK path below chains
              // zladiv's two); rounding differs from ZLARFG by an ulp
              const double rb = 1.0 / beta, dr = alpha.r - beta, den = dr * dr + alpha.i * alpha.i, rden = 1.0 / den;
              t = {(beta - alpha.r) * rb, -alpha.i * rb};
              scalx = {dr * rden, -alpha.i * rden};
            } else {
              t = {(beta - alpha.r) / beta, -alpha.i / beta};
              Z sc1 = zinv({alpha.r - beta, alpha.i});
              scalx = {sc1.r * skn, sc1.i * skn};
            }
            for (int k2 = 0; k2 < knt; ++k2) beta *= safmin;
          }
#pragma unroll
          for (int k = 0; k < KM; ++k) {
            const int r = (k * G + g) * SGS + lane;
            const Z vv = zmul(a[k], scalx);
            const Z vz = r < i ? Z{0.0, 0.0} : (r == i ? Z{1.0, 0.0} : vv);
            vb[vo + r] = sycl::double2(vz.r, vz.i);
            a[k] = r < i ? a[k] : (r == i ? Z{beta, 0.0} : vv);   // the final column (R part, beta, reflector), written after the barrier
          }
          if (g == 0 && lane == 0) { tau[2 * i] = t.r; tau[2 * i + 1] = t.i; tb[2 * buf] = t.r; tb[2 * buf + 1] = t.i; }
        }
#ifndef ABL_C
        it.barrier(sycl::access::fence_space::local_space);
#endif
        if (c == i) {
          double *dst = P + 2 * (int64_t)c * lda;
#pragma unroll
          for (int k = 0; k < KM; ++k) { const int r = (k * G + g) * SGS + lane; if (r < M) { dst[2 * r] = a[k].r; dst[2 * r + 1] = a[k].i; } }
        }
#ifdef QR2_TESTDELAY
        if (c == jb - 1) { double x_ = 1.0 + 1e-3 * (double)lane; for (int w_ = 0; w_ < pdly_; ++w_) x_ = x_ * 1.0000001 + 1e-9; if (x_ < 0.0) sx[0] = x_; }
#endif
        Z w = {0.0, 0.0};
        double tr = 0.0, ti = 0.0;
        bool gen = false;
        if (c > i) {
          tr = tb[2 * buf]; ti = tb[2 * buf + 1];
          gen = (tr != 0.0 || ti != 0.0);
          if (gen) {   // pass A: w = v^H a  (v kept in registers for pass B when it fits)
#pragma unroll
            for (int k = 0; k < KM; ++k) {
              const int r = (k * G + g) * SGS + lane;
              const sycl::double2 d = vb[vo + r];
              if (KEEPV) vreg[k] = {d.x(), d.y()};
              const Z p = zmulc(Z{d.x(), d.y()}, a[k]);
              w.r += p.r; w.i += p.i;
            }
#ifndef ABL_B
            red2(sg, w.r, w.i);
#endif
          }
        }
        if (G > 1) {   // the G sub-groups of a column combine their partial dot products
          if (gen && lane == 0) { wp[2 * sid] = w.r; wp[2 * sid + 1] = w.i; }
          it.barrier(sycl::access::fence_space::local_space);
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
              Z vk;
              if (KEEPV) vk = vreg[k];
              else { const int r = (k * G + g) * SGS + lane; const sycl::double2 d = vb[vo + r]; vk = {d.x(), d.y()}; }
              const Z p = zmul(f, vk);
              a[k].r -= p.r; a[k].i -= p.i;
            }
          }
          if (c == i + 1) publish(i + 1);   // next owner
        }
      }
    });
  });
}

// panel configuration: sub-groups per column G (1..4) and register slots KM (rows per lane). jb*G <= 32 sub-groups (512 work-items per group in 256-GRF mode).
// ALFDEVQR_G=1|2|3 forces G (experiments); default: see panel_g().
inline int panel_g_env() { static const int g = std::getenv("ALFDEVQR_G") ? std::atoi(std::getenv("ALFDEVQR_G")) : 0; return g; }
inline int panel_g(int64_t M, int jb) {
  int g = panel_g_env();
  if (g < 1) g = M <= 256 ? 1 : 2;
  while (g > 1 && jb * g > 32) --g;   // large-GRF mode: at most 32 sub-groups (512 work-items) per work-group
  return g;
}
inline int panel2_maxjb(int64_t M) { return M <= 768 ? 32 : (M <= 1536 ? 16 : 0); }
inline void panel2(sycl::queue &q, int M, int jb, double *P, int64_t lda, double *tau) {
  const int g = panel_g(M, jb);
#define PNL(GG, KMM) { panel2_t<GG, KMM>(q, M, jb, P, lda, tau); return; }
  if (g == 1) {
    if (M <= 64) PNL(1, 4) if (M <= 128) PNL(1, 8) if (M <= 192) PNL(1, 12) if (M <= 256) PNL(1, 16) if (M <= 384) PNL(1, 24) if (M <= 512) PNL(1, 32) if (M <= 640) PNL(1, 40) PNL(1, 48)
  } else if (g == 2) {
    if (M <= 128) PNL(2, 4) if (M <= 256) PNL(2, 8) if (M <= 384) PNL(2, 12) if (M <= 512) PNL(2, 16) if (M <= 640) PNL(2, 20) if (M <= 768) PNL(2, 24) if (M <= 1024) PNL(2, 32) if (M <= 1280) PNL(2, 40) PNL(2, 48)
  } else if (g == 3) {
    if (M <= 192) PNL(3, 4) if (M <= 384) PNL(3, 8) if (M <= 576) PNL(3, 12) if (M <= 768) PNL(3, 16) if (M <= 1024) PNL(3, 22) if (M <= 1280) PNL(3, 27) PNL(3, 32)
  } else {
    if (M <= 256) PNL(4, 4) if (M <= 512) PNL(4, 8) if (M <= 768) PNL(4, 12) PNL(4, 16)
  }
#undef PNL
}

// ---- apply: y <- (product of reflectors j in the given order) y  for every column y of C (TR: for every conj-transposed row of C) -----------
//   asc: j = 0..k-1, else k-1..0;  conjtau: use conj(tau_j).  Left,C = (asc, conj); Left,N = (desc, tau).  V: nq x k, unit diagonal implied, the
//   upper triangle is ignored. Row slots of sub-group g: (q*G+g)*16+lane. Columns nc (of y), NV per sub-group, G sub-groups per column group.
// Measured on the Max 1550 (phase 5): SLM bandwidth, not latency, bounds these kernels (a 16-byte-strided 8-byte load costs a full 256-byte pass of the SLM
// pipe), and so does global latency (~1 us per dependent load). Hence: the reflector is read from SLM as ONE 16-byte load per slot (double2) and reused for
// NV columns and for both passes (registers); the work-group stages CH reflectors (masked: 0 above the diagonal, 1 on it, 0 beyond nq) and their tau in
// SLM, double-buffered: the loads of chunk c+1 are issued into registers before chunk c is computed and stored to SLM after it; one barrier per chunk.
template <int Q, int NV, int G, bool TR, int CPW>
inline void apply2_t(sycl::queue &q, int nq, int k, int nc, bool asc, bool conjtau, const double *V, int64_t ldv, const double *tau, double *C, int64_t ldc) {
  constexpr int NQP = Q * G * SGS;                                 // padded rows
  constexpr int WI = CPW * G * SGS;
  constexpr int CHB = 65536 / (16 * 2 * NQP);
  constexpr int CH = CHB < 2 ? 2 : (CHB > 8 ? 8 : CHB);            // reflectors per chunk
  constexpr int NPF = (CH * NQP + WI - 1) / WI;                    // complex prefetch elements per work-item
  const int ncg = (nc + NV - 1) / NV;
  const int nwg = (ncg + CPW - 1) / CPW;
#ifdef QR2_TESTDELAY
  const int dly_ = G > 1 ? qr2_env_int("QR2_DELAY") : 0;
#endif
  q.submit([&](sycl::handler &h) {
    sycl::local_accessor<sycl::double2, 1> vs(sycl::range<1>(2 * (size_t)CH * NQP), h);   // [buf][u][row] complex
    sycl::local_accessor<double, 1> ts(sycl::range<1>(2 * 2 * CH), h);                    // [buf][u] tau
    sycl::local_accessor<double, 1> ex(sycl::range<1>(G > 1 ? 2 * 2 * (size_t)CPW * G * NV : 2), h);   // [parity][sub-group][v] partial dots
    sycl::ext::oneapi::experimental::properties props{sycl::ext::intel::experimental::grf_size<256>};
    h.parallel_for(sycl::nd_range<1>((size_t)nwg * WI, (size_t)WI), props, [=](sycl::nd_item<1> it) [[sycl::reqd_sub_group_size(SGS)]] {
      auto sg = it.get_sub_group();
      const int sid = (int)sg.get_group_id()[0];
      const int cg = sid / G, g = sid % G;
      const int lane = (int)sg.get_local_id()[0];
      const int tid = (int)it.get_local_linear_id();
      const int c0 = ((int)it.get_group(0) * CPW + cg) * NV;
      const int nch = (k + CH - 1) / CH;
      Z y[NV][Q];
#pragma unroll
      for (int v = 0; v < NV; ++v)
#pragma unroll
        for (int qq = 0; qq < Q; ++qq) {
          const int r = (qq * G + g) * SGS + lane;
          // branch-free: a load under a branch is waited for before the next one is issued; clamp the address and select
          const bool okc = (c0 + v < nc) && (r < nq);
          const int cc_ = okc ? c0 + v : 0, rr_ = okc ? r : 0;
          Z val;
          if (TR) { const double *p = C + 2 * ((int64_t)cc_ + (int64_t)rr_ * ldc); val = {p[0], -p[1]}; }
          else { const double *p = C + 2 * ((int64_t)rr_ + (int64_t)cc_ * ldc); val = {p[0], p[1]}; }
          if (!okc) val = {0.0, 0.0};
          y[v][qq] = val;
        }
      Z pf[NPF];
      Z tpf = {0.0, 0.0};
      auto prefetch = [&](int cc) {   // loads of chunk cc (steps cc*CH .. cc*CH+CH-1) into registers
#pragma unroll
        for (int p = 0; p < NPF; ++p) {
          const int e = tid + p * WI;
          const int u = e / NQP, r = e % NQP;
          const int s = cc * CH + u;
          const bool ok = (u < CH) && (s < k) && (r < nq);
          const int jr = asc ? s : k - 1 - s;
          const int j = ok ? jr : 0, rr = ok ? r : 0;
          const double *pv = V + 2 * ((int64_t)j * ldv + rr);
          Z val = {pv[0], pv[1]};   // always loaded (independent loads are pipelined)
          if (!ok || r < j) val = {0.0, 0.0}; else if (r == j) val = {1.0, 0.0};
          pf[p] = val;
        }
        const bool okt = tid < CH && cc * CH + tid < k;
        const int jt = okt ? (asc ? cc * CH + tid : k - 1 - (cc * CH + tid)) : 0;
        tpf = {tau[2 * jt], tau[2 * jt + 1]};
        if (!okt) tpf = {0.0, 0.0};
      };
      auto store = [&](int buf) {
        const size_t vo = (size_t)buf * CH * NQP;
#pragma unroll
        for (int p = 0; p < NPF; ++p) {
          const int e = tid + p * WI;
          if (e < CH * NQP) vs[vo + e] = sycl::double2(pf[p].r, pf[p].i);
        }
        if (tid < CH) { ts[2 * (buf * CH + tid)] = tpf.r; ts[2 * (buf * CH + tid) + 1] = tpf.i; }
      };
      prefetch(0);
      store(0);
      it.barrier(sycl::access::fence_space::local_space);
      int nex = 0;   // executed (non-skipped) steps: ex parity must alternate per EXECUTED step (review finding 1, 2026-11-12)
      for (int cc = 0; cc < nch; ++cc) {
        const int buf = cc & 1;
        if (cc + 1 < nch) prefetch(cc + 1);
        for (int u = 0; u < CH; ++u) {
          const int s = cc * CH + u;
          if (s >= k) break;
          Z t = {ts[2 * (buf * CH + u)], ts[2 * (buf * CH + u) + 1]};
          if (t.r == 0.0 && t.i == 0.0) continue;   // H = I (uniform: same tau for the whole work-group)
          if (conjtau) t.i = -t.i;
          const size_t vj = (size_t)buf * CH * NQP + (size_t)u * NQP;
          Z vr[Q];
          Z dot[NV];
#pragma unroll
          for (int v = 0; v < NV; ++v) dot[v] = {0.0, 0.0};
#pragma unroll
          for (int qq = 0; qq < Q; ++qq) {
            const sycl::double2 d = vs[vj + (qq * G + g) * SGS + lane];
            vr[qq] = {d.x(), d.y()};
#pragma unroll
            for (int v = 0; v < NV; ++v) { const Z p = zmulc(vr[qq], y[v][qq]); dot[v].r += p.r; dot[v].i += p.i; }
          }
#pragma unroll
          for (int v = 0; v < NV; ++v) red2(sg, dot[v].r, dot[v].i);
          if (G > 1) {
            const int eb = (0 * CPW * G + sid) * NV;
            if (lane == 0) {
#pragma unroll
              for (int v = 0; v < NV; ++v) { ex[2 * (eb + v)] = dot[v].r; ex[2 * (eb + v) + 1] = dot[v].i; }
            }
            it.barrier(sycl::access::fence_space::local_space);
#ifdef QR2_TESTDELAY
            if (sid == 1) { double x_ = 1.0 + 1e-3 * (double)lane; for (int w_ = 0; w_ < dly_; ++w_) x_ = x_ * 1.0000001 + 1e-9; if (x_ < 0.0) ex[0] = x_; }
#endif
#pragma unroll
            for (int v = 0; v < NV; ++v) {
              dot[v] = {0.0, 0.0};
#pragma unroll
              for (int gg = 0; gg < G; ++gg) { const int ei = (0 * CPW * G + cg * G + gg) * NV + v; dot[v].r += ex[2 * ei]; dot[v].i += ex[2 * ei + 1]; }
            }
          }
          Z coef[NV];
#pragma unroll
          for (int v = 0; v < NV; ++v) coef[v] = zmul(t, dot[v]);
#pragma unroll
          for (int qq = 0; qq < Q; ++qq) {
#pragma unroll
            for (int v = 0; v < NV; ++v) { const Z p = zmul(vr[qq], coef[v]); y[v][qq].r -= p.r; y[v][qq].i -= p.i; }
          }
          ++nex;
        }
        if (cc + 1 < nch) store(buf ^ 1);
        it.barrier(sycl::access::fence_space::local_space);
      }
#pragma unroll
      for (int v = 0; v < NV; ++v)
#pragma unroll
        for (int qq = 0; qq < Q; ++qq) {
          const int r = (qq * G + g) * SGS + lane;
          if (c0 + v < nc && r < nq) {
            if (TR) { double *p = C + 2 * ((int64_t)(c0 + v) + (int64_t)r * ldc); p[0] = y[v][qq].r; p[1] = -y[v][qq].i; }
            else { double *p = C + 2 * ((int64_t)r + (int64_t)(c0 + v) * ldc); p[0] = y[v][qq].r; p[1] = y[v][qq].i; }
          }
        }
    });
  });
}

constexpr int64_t APPLY_MAX_NQ = 1536;
// (Q row slots per lane, NV columns per sub-group, G sub-groups per column group) from nq: G*16*Q >= nq, NV*Q*4 GRF of the 256 hold y. Returns false if nq is too large.
template <bool TR>
inline bool apply_left(sycl::queue &q, int nq, int k, int nc, bool asc, bool conjtau, const double *V, int64_t ldv, const double *tau, double *C, int64_t ldc) {
  if (nc <= 0 || k <= 0) return true;
  if (nq > APPLY_MAX_NQ) return false;
#define AP(QQ, NVV, GG, CC) { apply2_t<QQ, NVV, GG, TR, CC>(q, nq, k, nc, asc, conjtau, V, ldv, tau, C, ldc); return true; }
  static const int anv = std::getenv("ALFDEVQR_ANV") ? std::atoi(std::getenv("ALFDEVQR_ANV")) : 2;
  if (nq <= 128) AP(8, 2, 1, 8)
  if (nq <= 288) AP(9, 2, 2, 4)
  if (nq <= 576) { if (anv == 4) AP(9, 4, 4, 4) else AP(9, 2, 4, 4) }
  if (nq <= 1152) { if (anv == 4) AP(9, 4, 8, 4) else AP(9, 2, 8, 4) }
  AP(12, 2, 8, 4)
#undef AP
}

static const int NB_SMALL = std::getenv("ALFDEVQR_NBS") ? std::atoi(std::getenv("ALFDEVQR_NBS")) : 32;
static const int NB_SMALL_M = std::getenv("ALFDEVQR_NBSM") ? std::atoi(std::getenv("ALFDEVQR_NBSM")) : 256;
// ---- drivers ---------------------------------------------------------------------------------------------------------
// ZGEQRF: panel2 + apply2 (no scratch). Returns 0, or -1 (not handled: m too large).
inline int geqrf2(sycl::queue &q, int64_t m, int64_t n, cd *A, int64_t lda, cd *tau, int nb) {
  const int64_t k = std::min(m, n);
  if (k == 0) return 0;
  if (nb < 1) nb = 1;
  if (nb > MAX_NB) nb = MAX_NB;
  if (m > MAX_M || panel2_maxjb(m) == 0) return -1;
  if (!q.is_in_order()) return -1;   // kernels below rely on in-order execution (no depends_on)
  {
    static const bool limits_ok = [&] {
      auto d = q.get_device();
      const size_t wg = d.get_info<sycl::info::device::max_work_group_size>(), slm = d.get_info<sycl::info::device::local_mem_size>();
      const auto sgs = d.get_info<sycl::info::device::sub_group_sizes>();
      const bool has16 = std::find(sgs.begin(), sgs.end(), (size_t)SGS) != sgs.end();
      const bool ok = wg >= 512 && slm >= 98304 && has16;
      if (!ok) std::fprintf(stderr, "alfdevqr: WARNING device limits too small for the device QR (max work-group %zu, SLM %zu B, sub-group 16: %d): using oneMKL\n", wg, slm, (int)has16);
      return ok;
    }();
    if (!limits_ok) return -1;
  }
  static const bool prof = std::getenv("ALFDEVQR_PROF") != nullptr;
  static double tp[2] = {0, 0}; static long ncall = 0;
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  double tl = prof ? now() : 0.0;
  auto lap = [&](int kk) { if (prof) { q.wait(); double t = now(); tp[kk] += t - tl; tl = t; } };
  for (int64_t j = 0, jbp; j < k; j += jbp) {
    const int Mj = (int)(m - j);
    const int nbj = (Mj <= NB_SMALL_M && nb == 16) ? NB_SMALL : nb;   // short panels: wider (fewer apply launches)
    const int jb = (int)std::min<int64_t>(std::min(nbj, panel2_maxjb(m - j)), k - j);
    jbp = jb;
    const int M = (int)(m - j);
    const int64_t nc = n - j - jb;
    cd *P = A + j + j * lda;
    panel2(q, M, jb, (double *)P, lda, (double *)(tau + j));
    lap(0);
    if (nc > 0) {
      if (!apply_left<false>(q, M, jb, (int)nc, true, true, (const double *)P, lda, (const double *)(tau + j), (double *)(A + j + (j + jb) * lda), lda)) return -2;
      lap(1);
    }
  }
  if (prof) { if (++ncall % 50 == 0) std::fprintf(stderr, "ALFDEVQR_PROF2 calls %ld  s total: panel %.3f apply %.3f  (m=%lld n=%lld nb=%d)\n", ncall, tp[0], tp[1], (long long)m, (long long)n, nb); }
  return 0;
}

// ZUNMQR. side 0 = Left, 1 = Right; trans 0 = N, else C (the LA_* convention of the shim: 0=N, 1=T, 2=C; T is not valid for complex).
// A holds k reflectors (nq x k, lda), nq = m (Left) or n (Right). 0 = done (asynchronous on the queue), -1 = declined before any change.
inline int unmqr2(sycl::queue &q, int side, int trans, int64_t m, int64_t n, int64_t k, const cd *A, int64_t lda, const cd *tau, cd *C, int64_t ldc) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;
  const int64_t nq = side == 0 ? m : n;
  if (nq > APPLY_MAX_NQ || k > nq) return -1;
  if (trans == 1) return -1;
  if (!q.is_in_order()) return -1;   // kernels below rely on in-order execution (no depends_on)
  if (side == 0) {   // Left,N: Q C = H_0 ... H_{k-1} C: descending, tau;  Left,C: Q^H C: ascending, conj(tau)
    const bool c = trans != 0;
    return apply_left<false>(q, (int)m, (int)k, (int)n, c, c, (const double *)A, lda, (const double *)tau, (double *)C, ldc) ? 0 : -1;
  }
  // Right: C Q = (Q^H C^H)^H, C Q^H = (Q C^H)^H: apply to the conj-transposed rows of C
  const bool tl = (trans == 0);   // Right,N -> Left,C (asc, conj tau);  Right,C -> Left,N (desc, tau)
  return apply_left<true>(q, (int)n, (int)k, (int)m, tl, tl, (const double *)A, lda, (const double *)tau, (double *)C, ldc) ? 0 : -1;
}
}  // namespace alfdevqr
