// ALF asynchronous submatrix-block path (GPU=intel, Phase 3b): the device sequence of one block of gpu_lowrank_mod
// (blk_begin gather, blk_end scale + flush gather + flush GEMMs, uploads) as SYCL kernels + oneMKL device GEMM, submitted
// on the sycl::queue of the OpenMP offload runtime (the one alf_la's shim uses; handed over by alfsh_get_queue), without
// waiting.  Only alfsb_wait blocks.  The queue must be in-order: that, and nothing else, orders the calls.
//
// Built as a SHARED library (libalfsb.so), because the device code of user SYCL kernels is only registered with the SYCL
// runtime when the sycl driver links it (an ifx-linked object gives "No kernel named ... was found").
// Every element-wise kernel computes the expressions of the OpenMP kernels of gpu_lowrank_mod (k_blk_gather, k_blk_scale,
// k_sub_flush_gather) in the same order; the GEMMs are the oneMKL calls of the dispatched la_gemm (same shapes/ops/scalars).
// All functions return 0, or -1 on a SYCL exception (printed on stderr).
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#ifndef ALFSB_PLANT_FUSE
#define ALFSB_PLANT_FUSE 0
#endif
using cd = std::complex<double>;
namespace {
sycl::queue *Q = nullptr;
sycl::event g_mark;  // completion event of the last alfsb_gather2 (waited for by alfsb_wait_mark)
// submission counters (ALFSB_COUNT=1: printed to stderr at exit; one increment per host submission through this library)
enum { C_WAIT, C_COPY, C_GATHER, C_SCALE, C_FGATHER, C_GEMM, C_SCALENS, C_FLUSH12, C_RANK2, C_GEMMB, C_CORR, C_GATHER2, C_MARK, C_FLUSH12R, C_SCALENS2, C_N };
const char *cname[C_N] = {"wait", "copy", "gather", "scale", "fgather", "gemm", "scale_ns", "flush12", "rank2", "gemm_batch", "corr", "gather2", "wait_mark", "flush12r", "scale_ns2"};
long cnt[C_N] = {0};
struct CounterReport {
  ~CounterReport() {
    const char *e = std::getenv("ALFSB_COUNT");
    if (!e || e[0] != '1') return;
    std::fprintf(stderr, "alfsb counts:");
    for (int i = 0; i < C_N; ++i) std::fprintf(stderr, " %s=%ld", cname[i], cnt[i]);
    std::fprintf(stderr, "\n");
  }
} counter_report;
template <class F> int guarded(const char *what, F &&f) {
  try {
    f();
    return 0;
  } catch (sycl::exception &e) {
    std::fprintf(stderr, "alfsb: %s: SYCL exception: %s\n", what, e.what());
  } catch (std::exception &e) {
    std::fprintf(stderr, "alfsb: %s: exception: %s\n", what, e.what());
  }
  return -1;
}
}  // namespace

extern "C" {
// q: sycl::queue* of the OpenMP runtime. Returns 0 if the queue is in-order, 1 if it is not (the caller must refuse), -1 error.
int alfsb_init(void *q) {
  if (!q) return -1;
  Q = static_cast<sycl::queue *>(q);
  return Q->is_in_order() ? 0 : 1;
}
void *alfsb_alloc_host(int64_t bytes) {
  try { return sycl::malloc_host(bytes, *Q); } catch (...) { return nullptr; }
}
void alfsb_free_host(void *p) { if (p && Q) sycl::free(p, *Q); }
int alfsb_wait() { ++cnt[C_WAIT]; return guarded("wait", [&] { Q->wait_and_throw(); }); }
int alfsb_copy(void *dst, const void *src, int64_t bytes) { ++cnt[C_COPY]; return guarded("copy", [&] { Q->memcpy(dst, src, (size_t)bytes); }); }

// gbk(a,c,fl) = Gk(bsk(a,fl), bsk(c,fl), nf), nf = nfa (fl=1) or nfb (fl=2): k_blk_gather
int alfsb_gather(int nfl, int nfa, int nfb, int ld, int b, int nd, cd *gbk, const int *bsk, const cd *Gk) {
  ++cnt[C_GATHER];
  return guarded("gather", [&] {
    Q->parallel_for(sycl::range<3>(nfl, b, b), [=](sycl::id<3> id) {
      long fl = id[0], c = id[1], a = id[2];
      long nf = (fl == 1) ? nfb : nfa;
      gbk[a + (long)nd * (c + (long)nd * fl)] =
          Gk[(bsk[a + (long)nd * fl] - 1) + (long)ld * ((bsk[c + (long)nd * fl] - 1) + (long)ld * (nf - 1))];
    });
  });
}
// k_blk_scale: rows, then columns of G(:,:,nf) scaled by the block's factors (cbk: za(1:b), zb(1:b) at nd+1:nd+b of each column)
int alfsb_scale(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd *Gk, const int *bsk, const cd *cbk) {
  ++cnt[C_SCALE];
  return guarded("scale", [&] {
    Q->parallel_for(sycl::range<3>(nfl, b, ld), [=](sycl::id<3> id) {
      long fl = id[0], a = id[1] + 1, j = id[2] + 1;
      long nf = (fl == 1) ? nfb : nfa;
      long sa = bsk[(a - 1) + (long)nd * fl];
      long cj = 0;
      for (long c = 1; c <= b; ++c)
        if (bsk[(c - 1) + (long)nd * fl] == j) cj = c;
      cd *Gf = Gk + (long)ld * ld * (nf - 1);
      cd v = cbk[(a - 1) + (long)ldc * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cbk[(nd + cj - 1) + (long)ldc * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cbk[(nd + a - 1) + (long)ldc * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
// k_sub_flush_gather: C(:,m) = G(:,S(m)); Qt(:,m) = e_S(m) - G(S(m),:)^T; slot(S(m),nf) = 0.  G: ld x ld of the flavour,
// slot: ld x nfg (column nf)
int alfsb_fgather(const cd *G, int ld, int nrow, const int *S, int k, int nf, cd *C, cd *Qt, int *slot) {
  ++cnt[C_FGATHER];
  return guarded("fgather", [&] {
    Q->parallel_for(sycl::range<2>(k, nrow), [=](sycl::id<2> id) {
      long m = id[0], i = id[1];
      long s = S[m];
      C[i + (long)ld * m] = G[i + (long)ld * (s - 1)];
      cd q = -G[(s - 1) + (long)ld * i];
      if (i + 1 == s) {
        q = q + 1.0;
        slot[i + (long)ld * (nf - 1)] = 0;
      }
      Qt[i + (long)ld * m] = q;
    });
  });
}
// zgemm, column major; ta/tb: 0 = N, 1 = T, 2 = C (the LA_* values)
int alfsb_gemm(int ta, int tb, int m, int n, int k, double ar, double ai, const cd *A, int lda, const cd *B, int ldb, double br,
               double bi, cd *C, int ldc) {
  auto t = [](int x) {
    return x == 0 ? oneapi::mkl::transpose::nontrans : (x == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans);
  };
  ++cnt[C_GEMM];
  return guarded("gemm", [&] {
    oneapi::mkl::blas::column_major::gemm(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, B, ldb, cd(br, bi), C, ldc);
  });
}


// ---- fused block steps (ALFSH_BLKFUSE, intel-batched step 2): fewer host submissions per update block ------------------------
// The block of gpu_lowrank_mod needs, per block, 13 submissions on the ASYNCBLK path.  The fused path packs the three host->device
// uploads of the block end into one copy from pinned memory and does the following in one kernel each:
//   alfsb_scale_ns : alfsb_scale + the nsigma_d copy (extra work-items at the end of the range); same per-element expressions
//   alfsb_flush12  : alfsb_fgather + the first flush GEMM (S2 = Qt M^T) of BOTH flavours; the GEMM is an own k x k product
//                    per row (sums over l in order), i.e. equal to oneMKL's result up to rounding, not bitwise
//   alfsb_rank2    : the second flush GEMM G -= X S2^T of both flavours (own register-blocked kernel, rounding-level differences)
// -DALFSB_PLANT_FUSE=1/2/3 plants one wrong index in flush12 / rank2 / scale_ns, 4/5 a relative 1e-9 error in flush12 / rank2, 6/7 a wrong index / 1e-9 error in gather2, 8/9 a missing row factor / 1e-9 error in flush12r (gate self-test only).
}  // extern "C"
namespace {
// complex multiply-accumulate with explicit FMAs (no __muldc3 NaN-recovery branches of std::complex operator*): (ar,ai) += (xr,xi) * (sr,si)
inline void cmac(double xr, double xi, double sr, double si, double &ar, double &ai) {
  ar = sycl::fma(xr, sr, ar);
  ar = sycl::fma(-xi, si, ar);
  ai = sycl::fma(xr, si, ai);
  ai = sycl::fma(xi, sr, ai);
}
}  // namespace
extern "C" {
int alfsb_scale_ns(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd *Gk, const int *bsk, const cd *cbk, int nsn,
                   const cd *nsrc, cd *ndst) {
  ++cnt[C_SCALENS];
  return guarded("scale_ns", [&] {
    const long nsc = (long)nfl * b * ld;
    Q->parallel_for(sycl::range<1>((size_t)(nsc + nsn)), [=](sycl::id<1> id) {
      long x = id[0];
      if (x >= nsc) {
#if ALFSB_PLANT_FUSE == 3
        ndst[(x - nsc + 1) % nsn] = nsrc[x - nsc];
#else
        ndst[x - nsc] = nsrc[x - nsc];
#endif
        return;
      }
      long j = x % ld + 1, a = (x / ld) % b + 1, fl = x / ((long)ld * b);
      long nf = (fl == 1) ? nfb : nfa;
      long sa = bsk[(a - 1) + (long)nd * fl];
      long cj = 0;
      for (long c = 1; c <= b; ++c)
        if (bsk[(c - 1) + (long)nd * fl] == j) cj = c;
      cd *Gf = Gk + (long)ld * ld * (nf - 1);
      cd v = cbk[(a - 1) + (long)ldc * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cbk[(nd + cj - 1) + (long)ldc * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cbk[(nd + a - 1) + (long)ldc * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
// Flush gathers + first GEMM of both flavours (fl = 0, 1 with nf = nf1, nf2, pending count kk1, kk2; kk = 0: flavour skipped).
// work-group (1, kmax, 16): m = local id 1 (column), i = 16 * tile + local id 2 (row, fastest: coalesced).
// Phase 1: item (m, i) with m < kk: q(i,m) = (e_S(m) - G(S(m),:)^T)(i) into local memory, C(i,m) = G(i,S(m)), slot reset.
// Phase 2: S2(i,m) = sum_l q(i,l) B(m,l), B = M packed at cbk + 2 nd + ldc fl (kk x kk, ld kk), the product S2 = Qt B^T of the GEMM.
// Columns kk <= m < kmax of C and S2 are set to zero (so that both flavours can be used with one common k in a batched GEMM).
int alfsb_flush12(int nfl, int nf1, int nf2, int kk1, int kk2, int kmax, int ld, int nrow, int nd, int ldc, const cd *G, int64_t sG,
                  cd *X, int64_t sX, const int *bS, const cd *cbk, cd *S2, int64_t sS2, int *slot) {
  ++cnt[C_FLUSH12];
  return guarded("flush12", [&] {
    const int ntile = (nrow + 15) / 16;
    Q->submit([&](sycl::handler &h) {
      sycl::local_accessor<cd, 1> qs(sycl::range<1>((size_t)kmax * 16), h);
      h.parallel_for(sycl::nd_range<3>(sycl::range<3>(nfl, kmax, (size_t)ntile * 16), sycl::range<3>(1, kmax, 16)),
                     [=](sycl::nd_item<3> it) {
        const int fl = it.get_group(0);
        const int kk = (fl == 0) ? kk1 : kk2;
        if (kk == 0) return;  // uniform over the work-group, before any barrier
        const int nf = (fl == 0) ? nf1 : nf2;
        const int m = it.get_local_id(1), il = it.get_local_id(2);
        const long i = (long)it.get_group(2) * 16 + il;
        const cd *Gf = G + sG * (nf - 1);
        cd *Xf = X + sX * (nf - 1);
        const int *S = bS + (long)nd * fl;
        const cd *B = cbk + 2 * nd + (long)ldc * fl;
        cd *S2f = S2 + sS2 * fl;
        if (m < kk && i < nrow) {
          long s = S[m];
          Xf[i + (long)ld * m] = Gf[i + (long)ld * (s - 1)];
          cd q = -Gf[(s - 1) + (long)ld * i];
          if (i + 1 == s) {
            q = q + 1.0;
            slot[i + (long)ld * (nf - 1)] = 0;
          }
          qs[m * 16 + il] = q;
        }
        it.barrier(sycl::access::fence_space::local_space);
        if (i < nrow) {
          if (m < kk) {
            double accr = 0, acci = 0;
            for (int l = 0; l < kk; ++l) {
#if ALFSB_PLANT_FUSE == 1
              const cd bb = B[l + (long)kk * m];
#else
              const cd bb = B[m + (long)kk * l];
#endif
              const cd qq = qs[l * 16 + il];
              cmac(qq.real(), qq.imag(), bb.real(), bb.imag(), accr, acci);
            }
#if ALFSB_PLANT_FUSE == 4
            accr *= 1.000000001; acci *= 1.000000001;  // subtle: relative 1e-9 error on S2
#endif
            S2f[i + (long)ld * m] = cd(accr, acci);
          } else {
            Xf[i + (long)ld * m] = 0;
            S2f[i + (long)ld * m] = 0;
          }
        }
      });
    });
  });
}
// G(nf) -= X(nf) S2^T for both flavours in one launch: G(i,j) = G(i,j) - sum_l X(i,l) S2(j,l), l < kk.
// work-group (1, 16, 16) = 16 columns x 16 rows, each item a 2 x 2 micro-tile (rows i, i + 16; columns j, j + 16); rows are the
// fastest index (coalesced). X and S2 come from L1/L2 (k <= 64, a few kB per work-group).
int alfsb_rank2(int nfl, int nf1, int nf2, int kk1, int kk2, int ld, cd *G, int64_t sG, const cd *X, int64_t sX, const cd *S2,
                int64_t sS2) {
  ++cnt[C_RANK2];
  return guarded("rank2", [&] {
    const int nt = (ld + 31) / 32;
    Q->parallel_for(sycl::nd_range<3>(sycl::range<3>(nfl, (size_t)nt * 16, (size_t)nt * 16), sycl::range<3>(1, 16, 16)),
                    [=](sycl::nd_item<3> it) {
      const int fl = it.get_group(0);
      const int kk = (fl == 0) ? kk1 : kk2;
      if (kk == 0) return;
      const int nf = (fl == 0) ? nf1 : nf2;
      const long i0 = (long)it.get_group(2) * 32 + it.get_local_id(2), i1 = i0 + 16;
      const long j0 = (long)it.get_group(1) * 32 + it.get_local_id(1), j1 = j0 + 16;
      const cd *Xf = X + sX * (nf - 1);
      const cd *Sf = S2 + sS2 * fl;
      double a00r = 0, a00i = 0, a01r = 0, a01i = 0, a10r = 0, a10i = 0, a11r = 0, a11i = 0;
      const bool vi0 = i0 < ld, vi1 = i1 < ld, vj0 = j0 < ld, vj1 = j1 < ld;
      for (int l = 0; l < kk; ++l) {
        cd x0 = vi0 ? Xf[i0 + (long)ld * l] : cd(0), x1 = vi1 ? Xf[i1 + (long)ld * l] : cd(0);
#if ALFSB_PLANT_FUSE == 2
        cd s0 = vj0 ? Sf[j0 + (long)ld * l] : cd(0), s1 = vj1 ? Sf[(j1 + 1) % ld + (long)ld * l] : cd(0);
#else
        cd s0 = vj0 ? Sf[j0 + (long)ld * l] : cd(0), s1 = vj1 ? Sf[j1 + (long)ld * l] : cd(0);
#endif
        cmac(x0.real(), x0.imag(), s0.real(), s0.imag(), a00r, a00i);
        cmac(x1.real(), x1.imag(), s0.real(), s0.imag(), a10r, a10i);
        cmac(x0.real(), x0.imag(), s1.real(), s1.imag(), a01r, a01i);
        cmac(x1.real(), x1.imag(), s1.real(), s1.imag(), a11r, a11i);
      }
      cd *Gf = G + sG * (nf - 1);
#if ALFSB_PLANT_FUSE == 5
      a00r *= 1.000000001; a00i *= 1.000000001;  // subtle: relative 1e-9 error on a quarter of the elements
#endif
      if (vi0 && vj0) Gf[i0 + (long)ld * j0] = cd(Gf[i0 + (long)ld * j0].real() - a00r, Gf[i0 + (long)ld * j0].imag() - a00i);
      if (vi1 && vj0) Gf[i1 + (long)ld * j0] = cd(Gf[i1 + (long)ld * j0].real() - a10r, Gf[i1 + (long)ld * j0].imag() - a10i);
      if (vi0 && vj1) Gf[i0 + (long)ld * j1] = cd(Gf[i0 + (long)ld * j1].real() - a01r, Gf[i0 + (long)ld * j1].imag() - a01i);
      if (vi1 && vj1) Gf[i1 + (long)ld * j1] = cd(Gf[i1 + (long)ld * j1].real() - a11r, Gf[i1 + (long)ld * j1].imag() - a11i);
    });
  });
}

// Level 5 (ALFSH_BLKFUSE=5): the flush gathers read the UNSCALED G and apply the block's scale factors on the fly (inv: site -> index in the block,
// 0 if not a block site; per flavour, ld ints), so the scale kernel is no longer on the path to the flush factors and the prefetch kernel:
// order flush12r, gather2, scale_ns2, GEMM 2.  Same operations in the same order as scale + fgather (std::complex products, bitwise equal X and Qt).
// alfsb_scale_ns2 = alfsb_scale_ns with the search over the block sites replaced by the table lookup (cj = inv[j]).
int alfsb_scale_ns2(int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd *Gk, const int *bsk, const cd *cbk, int nsn,
                    const cd *nsrc, cd *ndst, const int *inv) {
  ++cnt[C_SCALENS2];
  return guarded("scale_ns2", [&] {
    const long nsc = (long)nfl * b * ld;
    Q->parallel_for(sycl::range<1>((size_t)(nsc + nsn)), [=](sycl::id<1> id) {
      long x = id[0];
      if (x >= nsc) {
        ndst[x - nsc] = nsrc[x - nsc];
        return;
      }
      long j = x % ld + 1, a = (x / ld) % b + 1, fl = x / ((long)ld * b);
      long nf = (fl == 1) ? nfb : nfa;
      long sa = bsk[(a - 1) + (long)nd * fl];
      long cj = inv[(j - 1) + (long)ld * fl];
      cd *Gf = Gk + (long)ld * ld * (nf - 1);
      cd v = cbk[(a - 1) + (long)ldc * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cbk[(nd + cj - 1) + (long)ldc * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cbk[(nd + a - 1) + (long)ldc * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
int alfsb_flush12r(int nfl, int nf1, int nf2, int kk1, int kk2, int kmax, int ld, int nrow, int nd, int ldc, const cd *G, int64_t sG,
                   cd *X, int64_t sX, const int *bS, const cd *cbk, cd *S2, int64_t sS2, int *slot, const int *inv) {
  ++cnt[C_FLUSH12R];
  return guarded("flush12r", [&] {
    const int ntile = (nrow + 15) / 16;
    Q->submit([&](sycl::handler &h) {
      sycl::local_accessor<cd, 1> qs(sycl::range<1>((size_t)kmax * 16), h);
      h.parallel_for(sycl::nd_range<3>(sycl::range<3>(nfl, kmax, (size_t)ntile * 16), sycl::range<3>(1, kmax, 16)),
                     [=](sycl::nd_item<3> it) {
        const int fl = it.get_group(0);
        const int kk = (fl == 0) ? kk1 : kk2;
        if (kk == 0) return;
        const int nf = (fl == 0) ? nf1 : nf2;
        const int m = it.get_local_id(1), il = it.get_local_id(2);
        const long i = (long)it.get_group(2) * 16 + il;
        const cd *Gf = G + sG * (nf - 1);
        cd *Xf = X + sX * (nf - 1);
        const int *S = bS + (long)nd * fl;
        const int *iv = inv + (long)ld * fl;
        const cd *cb = cbk + (long)ldc * fl;
        const cd *B = cbk + 2 * nd + (long)ldc * fl;
        cd *S2f = S2 + sS2 * fl;
        if (m < kk && i < nrow) {
          const long s = S[m];
          const long am = iv[s - 1];  // block index (1-based) of the pending site, > 0
          const long ai = iv[i];      // block index of row i, 0 if not a block site
          // C(i,m) = G'(i,s): row factor of i (if i in the block), then the column factor of s
          cd g = Gf[i + (long)ld * (s - 1)];
#if ALFSB_PLANT_FUSE != 8
          if (ai > 0) g = cb[ai - 1] * g;
#endif
          g = cb[nd + am - 1] * g;
          Xf[i + (long)ld * m] = g;
          // q = e_s - G'(s,:)^T: row factor of s, then (if i in the block) the column factor of i
          cd v = cb[am - 1] * Gf[(s - 1) + (long)ld * i];
          if (ai > 0) v = cb[nd + ai - 1] * v;
          cd q = -v;
          if (i + 1 == s) {
            q = q + 1.0;
            slot[i + (long)ld * (nf - 1)] = 0;
          }
          qs[m * 16 + il] = q;
        }
        it.barrier(sycl::access::fence_space::local_space);
        if (i < nrow) {
          if (m < kk) {
            double accr = 0, acci = 0;
            for (int l = 0; l < kk; ++l) {
              const cd bb = B[m + (long)kk * l];
              const cd qq = qs[l * 16 + il];
              cmac(qq.real(), qq.imag(), bb.real(), bb.imag(), accr, acci);
            }
#if ALFSB_PLANT_FUSE == 9
            accr *= 1.000000001; acci *= 1.000000001;  // subtle: relative 1e-9 error on S2
#endif
            S2f[i + (long)ld * m] = cd(accr, acci);
          } else {
            Xf[i + (long)ld * m] = 0;
            S2f[i + (long)ld * m] = 0;
          }
        }
      });
    });
  });
}
// Prefetch of the next block's principal block (ALFSH_BLKPF=1): G'(B',B') = G(B',B') - X(B',:) S2(B',:)^T of the flush that is still to be
// done (sites B' disjoint from the block sites, so the scaling does not touch them), written straight into host memory, and the next block's
// site list bsn copied into bsd for its own scale.  The completion event is kept: alfsb_wait_mark waits for this kernel only, not for the
// flush GEMMs queued after it.  The products are the sums of alfsb_rank2 (same order, FMAs, G - acc).
int alfsb_gather2(int nfl, int nf1, int nf2, int kk1, int kk2, int ld, int b, int nd, cd *gbk, int *bsd, const int *bsn, const cd *G,
                  int64_t sG, const cd *X, int64_t sX, const cd *S2, int64_t sS2) {
  ++cnt[C_GATHER2];
  return guarded("gather2", [&] {
    g_mark = Q->parallel_for(sycl::range<3>(nfl, b, b), [=](sycl::id<3> id) {
      const long fl = id[0], c = id[1], a = id[2];
      const int nf = (fl == 0) ? nf1 : nf2;
      const int kk = (fl == 0) ? kk1 : kk2;
#if ALFSB_PLANT_FUSE == 6
      const long i = bsn[a + (long)nd * fl] - 1, j = bsn[(c + 1) % b + (long)nd * fl] - 1;
#else
      const long i = bsn[a + (long)nd * fl] - 1, j = bsn[c + (long)nd * fl] - 1;
#endif
      const cd *Gf = G + sG * (nf - 1);
      double vr = Gf[i + (long)ld * j].real(), vi = Gf[i + (long)ld * j].imag();
      if (kk > 0) {
        const cd *Xf = X + sX * (nf - 1);
        const cd *Sf = S2 + sS2 * fl;
        double ar = 0, ai = 0;
        for (int l = 0; l < kk; ++l) {
          const cd x = Xf[i + (long)ld * l], sv = Sf[j + (long)ld * l];
          cmac(x.real(), x.imag(), sv.real(), sv.imag(), ar, ai);
        }
        vr -= ar;
        vi -= ai;
      }
#if ALFSB_PLANT_FUSE == 7
      vr *= 1.000000001; vi *= 1.000000001;  // subtle: relative 1e-9 error on the prefetched block
#endif
      gbk[a + (long)nd * (c + (long)nd * fl)] = cd(vr, vi);
      if (a == 0 && c == 0)
        for (int q = 0; q < nd; ++q) bsd[q + (long)nd * fl] = bsn[q + (long)nd * fl];
    });
  });
}
int alfsb_wait_mark() {
  ++cnt[C_MARK];
  return guarded("wait_mark", [&] { g_mark.wait_and_throw(); });
}

// ---- observables: SYCL port of the HIP k_corr (alf_hip_kernels.hip.cpp), switch ALFSH_OBSCORR=1 in gpu_obs_mod ----------------
// (Prog/gpu_obs_mod.F90 corr; not called unless the switch is set. -DALFSB_PLANT_CORR plants two wrong indices (kind 2 and kind 4): gate self-test only.)
// One work-group (256) per output m = (imj,no,no2); work-groups MT..MT+Norb-1 do the acc0 sums (kinds 3, 4). Same expressions as the
// OpenMP kernels of gpu_obs_mod corr; the order of the sum over I differs (work-group tree reduction), i.e. rounding-level changes.
// Callable without the block path: q is the OpenMP interop queue (alfsh_get_queue); waits at the end (acc is read by OpenMP code).
int alfsb_corr(void *q, int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld, const cd *G,
               const int *jc, const int *inv, cd *acc, cd *acc0) {
  if (!Q) {
    if (!q) return -1;
    Q = static_cast<sycl::queue *>(q);
  }
  ++cnt[C_CORR];
  return guarded("corr", [&] {
    const long MT = (long)N * Norb * Norb;
    const long nwg = MT + ((kind == 3 || kind == 4) ? Norb : 0);
    constexpr int RB = 256;
    Q->parallel_for(sycl::nd_range<1>((size_t)nwg * RB, RB), [=](sycl::nd_item<1> it) {
      const long m = it.get_group(0);
      const bool zero_part = (m >= MT);
      int imj = 0, no, no2 = 0;
      if (!zero_part) { imj = m % N; no = (m / N) % Norb; no2 = m / ((long)N * Norb); }
      else no = (int)(m - MT);
      cd s1 = 0, s2 = 0, s3 = 0;
      const long ld2 = (long)ld * ld;
      const cd one(1.0, 0.0);
      for (int I = it.get_local_id(0); I < N; I += RB) {
        if (!zero_part) {
          int J = jc[I + (long)imj * N] - 1;
          int i1 = inv[I + (long)no * N] - 1, j1 = inv[J + (long)no2 * N] - 1;
          if (kind == 1) {
            cd z = 0;
            for (int nf = 0; nf < NFL; nf++) { cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one; z = z + a; }
            s1 = s1 + z * ns;
          } else if (kind == 2) {
            cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
#ifdef ALFSB_PLANT_CORR
            s1 = s1 + (a * G[j1 + (long)i1 * ld]) * ns;
#else
            s1 = s1 + (a * G[i1 + (long)j1 * ld]) * ns;
#endif
          } else if (kind == 3) {
            cd zi = 0, zj = 0, z = 0;
            for (int nf = 0; nf < NFL; nf++) {
              zi = zi + (one - G[i1 + (long)i1 * ld + nf * ld2]);
              zj = zj + (one - G[j1 + (long)j1 * ld + nf * ld2]);
              cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one;
              z = z + a * G[i1 + (long)j1 * ld + nf * ld2];
            }
            s1 = s1 + (zi * ns) * (zj * ns);
            s1 = s1 + z * ns;
          } else {
            cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
            cd b = -G[j1 + (long)i1 * ld + ld2]; if (i1 == j1) b = b + one;
            cd zxy = a * G[i1 + (long)j1 * ld + ld2] + b * G[i1 + (long)j1 * ld];
            cd zi = (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
            cd zj = (one - G[j1 + (long)j1 * ld + ld2]) - (one - G[j1 + (long)j1 * ld]);
#ifdef ALFSB_PLANT_CORR
            cd zz = a * G[i1 + (long)j1 * ld] + b * G[j1 + (long)i1 * ld + ld2] + zi * zj;
#else
            cd zz = a * G[i1 + (long)j1 * ld] + b * G[i1 + (long)j1 * ld + ld2] + zi * zj;
#endif
            s1 = s1 + zz; s2 = s2 + zxy; s3 = s3 + (2.0 * zxy + zz) / 3.0;
          }
        } else {
          int i1 = inv[I + (long)no * N] - 1;
          if (kind == 3) {
            cd zi = 0;
            for (int nf = 0; nf < NFL; nf++) zi = zi + (one - G[i1 + (long)i1 * ld + nf * ld2]);
            s1 = s1 + zi * ns;
          } else {
            s1 = s1 + ((one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]));
          }
        }
      }
      auto g = it.get_group();
      double r0 = sycl::reduce_over_group(g, s1.real(), sycl::plus<double>());
      double r1 = sycl::reduce_over_group(g, s1.imag(), sycl::plus<double>());
      double r2 = sycl::reduce_over_group(g, s2.real(), sycl::plus<double>());
      double r3 = sycl::reduce_over_group(g, s2.imag(), sycl::plus<double>());
      double r4 = sycl::reduce_over_group(g, s3.real(), sycl::plus<double>());
      double r5 = sycl::reduce_over_group(g, s3.imag(), sycl::plus<double>());
      if (it.get_local_id(0) == 0) {
        const cd w(wr, wi);
        const cd t1(r0, r1), t2(r2, r3), t3(r4, r5);
        if (!zero_part) {
          acc[m + (long)(k - 1) * MT] += t1 * w;
          if (kind == 4) { acc[m + (long)k * MT] += t2 * w; acc[m + (long)(k + 1) * MT] += t3 * w; }
        } else {
          acc0[no + (long)(k - 1) * Norb] += t1 * w;
          if (kind == 4) acc0[no + (long)(k + 1) * Norb] += (t1 * w) / sq3;
        }
      }
    });
    Q->wait_and_throw();
  });
}


// ---- multi-walker (batched) variants, intel-batched branch, NOT yet called by Fortran (docs/gpu/platforms/INTEL_BATCHED.md) ----
// W walkers, walker w's arrays at base + w * stride (strides in elements). Per element, every kernel computes exactly the
// expression of the single-walker kernel above (same operands, same order), so walker w of a batched launch is bitwise
// equal to the single-walker launch on walker w's data. Per-walker scalars (nfa/nfb, k, nf) are the same for all walkers in
// lockstep (same model, same slice and block); a walker with nothing to do (k = 0 in the flush) is skipped via kw[w] = 0.
int alfsb_gather_w(int W, int nfl, int nfa, int nfb, int ld, int b, int nd, cd *gbk, int64_t sgb, const int *bsk, int64_t sbs,
                   const cd *Gk, int64_t sG) {
  return guarded("gather_w", [&] {
    Q->parallel_for(sycl::range<3>((size_t)W * nfl, b, b), [=](sycl::id<3> id) {
      long w = id[0] / nfl, fl = id[0] % nfl, c = id[1], a = id[2];
      long nf = (fl == 1) ? nfb : nfa;
      const int *bs = bsk + sbs * w;
      const cd *G = Gk + sG * w;
      gbk[sgb * w + a + (long)nd * (c + (long)nd * fl)] =
          G[(bs[a + (long)nd * fl] - 1) + (long)ld * ((bs[c + (long)nd * fl] - 1) + (long)ld * (nf - 1))];
    });
  });
}
int alfsb_scale_w(int W, int nfl, int nfa, int nfb, int ld, int b, int nd, int ldc, cd *Gk, int64_t sG, const int *bsk, int64_t sbs,
                  const cd *cbk, int64_t scb) {
  return guarded("scale_w", [&] {
    Q->parallel_for(sycl::range<3>((size_t)W * nfl, b, ld), [=](sycl::id<3> id) {
      long w = id[0] / nfl, fl = id[0] % nfl, a = id[1] + 1, j = id[2] + 1;
      long nf = (fl == 1) ? nfb : nfa;
      const int *bs = bsk + sbs * w;
      const cd *cb = cbk + scb * w;
      long sa = bs[(a - 1) + (long)nd * fl];
      long cj = 0;
      for (long c = 1; c <= b; ++c)
        if (bs[(c - 1) + (long)nd * fl] == j) cj = c;
      cd *Gf = Gk + sG * w + (long)ld * ld * (nf - 1);
      cd v = cb[(a - 1) + (long)ldc * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cb[(nd + cj - 1) + (long)ldc * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cb[(nd + a - 1) + (long)ldc * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
// kmax = max over walkers of k; kw[w] = k of walker w (device array), columns m >= kw[w] are not touched
int alfsb_fgather_w(int W, const cd *G, int64_t sG, int ld, int nrow, const int *S, int64_t sS, const int *kw, int kmax, int nf,
                    cd *C, int64_t sC, cd *Qt, int64_t sQ, int *slot, int64_t ssl) {
  return guarded("fgather_w", [&] {
    Q->parallel_for(sycl::range<2>((size_t)W * kmax, nrow), [=](sycl::id<2> id) {
      long w = id[0] / kmax, m = id[0] % kmax, i = id[1];
      if (m >= kw[w]) return;
      const cd *Gw = G + sG * w;
      long s = S[sS * w + m];
      C[sC * w + i + (long)ld * m] = Gw[i + (long)ld * (s - 1)];
      cd q = -Gw[(s - 1) + (long)ld * i];
      if (i + 1 == s) {
        q = q + 1.0;
        slot[ssl * w + i + (long)ld * (nf - 1)] = 0;
      }
      Qt[sQ * w + i + (long)ld * m] = q;
    });
  });
}
// W gemms with walker strides. mode 0: W separate oneMKL gemm calls (bitwise equal to the single-walker path, one host
// submission each); mode 1: one oneMKL gemm_batch (strided) call (one submission; NOT bitwise equal to the single call for
// most N, Phase-6 job 560787). Same k for all walkers (lockstep: walkers with smaller k are padded by the caller).
int alfsb_gemm_w(int mode, int W, int ta, int tb, int m, int n, int k, double ar, double ai, const cd *A, int lda, int64_t sA,
                 const cd *B, int ldb, int64_t sB, double br, double bi, cd *C, int ldc, int64_t sC) {
  auto t = [](int x) {
    return x == 0 ? oneapi::mkl::transpose::nontrans : (x == 1 ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::conjtrans);
  };
  return guarded("gemm_w", [&] {
    if (mode == 0)
      for (int w = 0; w < W; ++w)
        oneapi::mkl::blas::column_major::gemm(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A + sA * w, lda, B + sB * w, ldb, cd(br, bi),
                                              C + sC * w, ldc);
    else {
      ++cnt[C_GEMMB];
      oneapi::mkl::blas::column_major::gemm_batch(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, sA, B, ldb, sB, cd(br, bi), C, ldc,
                                                  sC, W);
    }
  });
}
}
