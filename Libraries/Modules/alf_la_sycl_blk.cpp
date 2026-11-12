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
using cd = std::complex<double>;
namespace {
sycl::queue *Q = nullptr;
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
}
