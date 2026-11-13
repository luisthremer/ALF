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
#if defined(__has_include)
#if __has_include(<sycl/ext/oneapi/experimental/enqueue_functions.hpp>)
#include <sycl/ext/oneapi/experimental/enqueue_functions.hpp>
#define ALFSB_HAVE_ENQ 1
#endif
#endif
using cd = std::complex<double>;
namespace {
sycl::queue *Q = nullptr;
// BLKPACK path (ALFSH_BLKPACK=1 in gpu_lowrank_mod): kernels / copy are submitted with the enqueue functions (no sycl::event is created:
// cheaper host submission).  ALFSB_NOENQ=1 falls back to queue::parallel_for / memcpy (ablation).
int g_enq = -1;
bool use_enq() {
#ifdef ALFSB_HAVE_ENQ
  if (g_enq < 0) { const char *e = std::getenv("ALFSB_NOENQ"); g_enq = (e && e[0] == '1') ? 0 : 1; }
  return g_enq == 1;
#else
  return false;
#endif
}
template <int D, class K> void launch(sycl::range<D> r, K k) {
#ifdef ALFSB_HAVE_ENQ
  if (use_enq()) { sycl::ext::oneapi::experimental::parallel_for(*Q, r, k); return; }
#endif
  Q->parallel_for(r, k);
}
constexpr int MAXND = 64;          // block size capacity of the by-value block-site argument of the packed gather
struct IntBuf { int v[2 * MAXND]; };
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
void *alfsb_alloc_dev(int64_t bytes) {
  try { return sycl::malloc_device(bytes, *Q); } catch (...) { return nullptr; }
}
void alfsb_free_dev(void *p) { if (p && Q) sycl::free(p, *Q); }
int alfsb_wait() { return guarded("wait", [&] { Q->wait_and_throw(); }); }
int alfsb_copy(void *dst, const void *src, int64_t bytes) { return guarded("copy", [&] { Q->memcpy(dst, src, (size_t)bytes); }); }

// gbk(a,c,fl) = Gk(bsk(a,fl), bsk(c,fl), nf), nf = nfa (fl=1) or nfb (fl=2): k_blk_gather
int alfsb_gather(int nfl, int nfa, int nfb, int ld, int b, int nd, cd *gbk, const int *bsk, const cd *Gk) {
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
  return guarded("gemm", [&] {
    oneapi::mkl::blas::column_major::gemm(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, B, ldb, cd(br, bi), C, ldc);
  });
}

// ---- BLKPACK variants (bitwise equal to the functions above; fewer / cheaper submissions) ---------------------------------------
// host -> device copy of the packed block data (one copy instead of four), no event
int alfsb_upload(void *dst, const void *src, int64_t bytes) {
  return guarded("upload", [&] {
#ifdef ALFSB_HAVE_ENQ
    if (use_enq()) { sycl::ext::oneapi::experimental::memcpy(*Q, dst, src, (size_t)bytes); return; }
#endif
    Q->memcpy(dst, src, (size_t)bytes);
  });
}
// alfsb_gather with (a) the block sites bs_h(nd,2) (HOST array, ints) passed BY VALUE (no upload of bsd) and (b) the result written straight
// into gbk, which may be host USM (no device-to-host copy afterwards).  Same element expression.  Returns -2 if nd > MAXND.
int alfsb_gatherv(int nfl, int nfa, int nfb, int ld, int b, int nd, cd *gbk, const int *bs_h, const cd *Gk) {
  if (nd > MAXND) return -2;
  IntBuf bsv;
  for (int i = 0; i < 2 * MAXND; ++i) bsv.v[i] = (i < 2 * nd) ? bs_h[i] : 0;
  static int gdev = -1;                 // ablation ALFSB_GDEV=1: the kernel writes a device buffer, then one copy to gbk (as the unpacked path does)
  static cd *gbuf = nullptr;
  if (gdev < 0) { const char *e = std::getenv("ALFSB_GDEV"); gdev = (e && e[0] == '1') ? 1 : 0; }
  cd *gout = gbk;
  if (gdev == 1) {
    if (!gbuf) gbuf = (cd *)sycl::malloc_device((size_t)2 * nd * nd * sizeof(cd), *Q);
    gout = gbuf;
  }
  int rc = guarded("gatherv", [&] {
    launch<3>(sycl::range<3>(nfl, b, b), [=](sycl::item<3> it) {
      long fl = it.get_id(0), c = it.get_id(1), a = it.get_id(2);
      long nf = (fl == 1) ? nfb : nfa;
      gout[a + (long)nd * (c + (long)nd * fl)] =
          Gk[(bsv.v[a + nd * fl] - 1) + (long)ld * ((bsv.v[c + nd * fl] - 1) + (long)ld * (nf - 1))];
    });
  });
  if (rc == 0 && gdev == 1) rc = guarded("gatherv copy", [&] { Q->memcpy(gbk, gbuf, (size_t)2 * nd * nd * sizeof(cd)); });
  return rc;
}
// alfsb_scale on the packed data (bsk: b ints per flavour; cbk: za(1:b), zb(1:b) per flavour, leading dimension ldz = 2 b); the nns field values
// nss -> nsd (nsigma_d(vmin:vmax,t)) are copied by the work items (fl = 0, a = 1, j <= nns) of the same launch (the scale kernel does not read nsigma_d;
// nns <= b <= ld).  Same 3-D launch and per-element expressions as alfsb_scale.
int alfsb_scalep(int nfl, int nfa, int nfb, int ld, int b, int ldz, cd *Gk, const int *bsk, const cd *cbk, cd *nsd, const cd *nss, int nns) {
  return guarded("scalep", [&] {
    launch<3>(sycl::range<3>(nfl, b, ld), [=](sycl::item<3> it) {
      long fl = it.get_id(0), a = it.get_id(1) + 1, j = it.get_id(2) + 1;
      if (fl == 0 && a == 1 && j <= nns) nsd[j - 1] = nss[j - 1];
      long nf = (fl == 1) ? nfb : nfa;
      long sa = bsk[(a - 1) + (long)b * fl];
      long cj = 0;
      for (long c = 1; c <= b; ++c)
        if (bsk[(c - 1) + (long)b * fl] == j) cj = c;
      cd *Gf = Gk + (long)ld * ld * (nf - 1);
      cd v = cbk[(a - 1) + (long)ldz * fl] * Gf[(sa - 1) + (long)ld * (j - 1)];
      if (cj > 0) v = cbk[(b + cj - 1) + (long)ldz * fl] * v;
      Gf[(sa - 1) + (long)ld * (j - 1)] = v;
      if (cj == 0) Gf[(j - 1) + (long)ld * (sa - 1)] = cbk[(b + a - 1) + (long)ldz * fl] * Gf[(j - 1) + (long)ld * (sa - 1)];
    });
  });
}
// alfsb_fgather of both flavours in ONE launch (nfl = 1 or 2; flavour fl: G_fl (ld x ld), S_fl (kk_fl ints, device), C_fl, Qt_fl, slot column nf_fl).
// Each flavour has its own Qt buffer here (the sequential version reused one), the values are the same.  kk_fl = 0: flavour skipped.
int alfsb_fgather2(int nfl, int ld, int nrow, int kk0, int kk1, int nf0, int nf1, const cd *G0, const cd *G1, const int *S0, const int *S1,
                   cd *C0, cd *C1, cd *Q0, cd *Q1, int *slot) {
  int kmax = kk0 > kk1 ? kk0 : kk1;
  if (nfl < 2) kmax = kk0;
  if (kmax <= 0) return 0;
  return guarded("fgather2", [&] {
    launch<3>(sycl::range<3>(nfl, kmax, nrow), [=](sycl::item<3> it) {
      long fl = it.get_id(0), m = it.get_id(1), i = it.get_id(2);
      int kk = fl ? kk1 : kk0;
      if (m >= kk) return;
      const cd *G = fl ? G1 : G0;
      const int *S = fl ? S1 : S0;
      cd *C = fl ? C1 : C0;
      cd *Qt = fl ? Q1 : Q0;
      int nf = fl ? nf1 : nf0;
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
// ---- BLKPACK=2 (NOT bitwise equal: the k-sum of Y = Qt M^T is accumulated in the order j = 0..kk-1 by this kernel, not by oneMKL's zgemm) -----
// alfsb_fgather2 + the first GEMM of the flush (Y = Qt M^T, Qt = E_S - G(S,:)^T, M kk x kk with leading dimension kk) in ONE launch; Qt is not stored.
// Work item (fl, i, m): C(i,m) = G(i,S(m)); slot reset; Y(i,m) = sum_j Qt(i,j) M(m,j).  kk_fl = 0: flavour skipped.
int alfsb_fgemm1(int nfl, int ld, int nrow, int kk0, int kk1, int nf0, int nf1, const cd *G0, const cd *G1, const int *S0, const int *S1,
                 const cd *M0, const cd *M1, cd *C0, cd *C1, cd *Y0, cd *Y1, int *slot) {
  int kmax = kk0 > kk1 ? kk0 : kk1;
  if (nfl < 2) kmax = kk0;
  if (kmax <= 0) return 0;
  return guarded("fgemm1", [&] {
    launch<3>(sycl::range<3>(nfl, nrow, kmax), [=](sycl::item<3> it) {
      long fl = it.get_id(0), i = it.get_id(1), m = it.get_id(2);
      int kk = fl ? kk1 : kk0;
      if (m >= kk) return;
      const cd *G = fl ? G1 : G0;
      const int *S = fl ? S1 : S0;
      const cd *M = fl ? M1 : M0;
      cd *C = fl ? C1 : C0;
      cd *Y = fl ? Y1 : Y0;
      int nf = fl ? nf1 : nf0;
      long s = S[m];
      C[i + (long)ld * m] = G[i + (long)ld * (s - 1)];
      if (i + 1 == s) slot[i + (long)ld * (nf - 1)] = 0;
      cd acc(0.0, 0.0);
      for (int j = 0; j < kk; ++j) {
        long sj = S[j];
        cd q = -G[(sj - 1) + (long)ld * i];
        if (i + 1 == sj) q = q + 1.0;
        acc += q * M[m + (long)kk * j];
      }
      Y[i + (long)ld * m] = acc;
    });
  });
}
}
