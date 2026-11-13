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

#ifdef ALFSB_WIP_CORR  // WORK IN PROGRESS: not built, not tested (docs/gpu/platforms/INTEL_BATCHED.md section 5); compiled only with -DALFSB_WIP_CORR
// ---- observables: SYCL port of the HIP k_corr (alf_hip_kernels.hip.cpp), switch ALFSH_OBSCORR=1 in gpu_obs_mod ----------------
// One work-group (128) per output m = (imj,no,no2); work-groups MT..MT+Norb-1 do the acc0 sums (kinds 3, 4). Same expressions as the
// OpenMP kernels of gpu_obs_mod corr; the order of the sum over I differs (work-group tree reduction), i.e. rounding-level changes.
// Callable without the block path: q is the OpenMP interop queue (alfsh_get_queue); waits at the end (acc is read by OpenMP code).
int alfsb_corr(void *q, int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld, const cd *G,
               const int *jc, const int *inv, cd *acc, cd *acc0) {
  if (!Q) {
    if (!q) return -1;
    Q = static_cast<sycl::queue *>(q);
  }
  return guarded("corr", [&] {
    const long MT = (long)N * Norb * Norb;
    const long nwg = MT + ((kind == 3 || kind == 4) ? Norb : 0);
    constexpr int RB = 128;
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
            s1 = s1 + (a * G[i1 + (long)j1 * ld]) * ns;
          } else if (kind == 3) {
            cd zi = 0, zj = 0, z = 0;
            for (int nf = 0; nf < NFL; nf++) {
              zi = zi + (one - G[i1 + (long)i1 * ld + nf * ld2]);
              zj = zj + (one - G[j1 + (long)j1 * ld + nf * ld2]);
              cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one;
              z = z + a * G[i1 + (long)j1 * ld + nf * ld2];
            }
            s1 = s1 + ((zi * ns) * (zj * ns) + z * ns);
          } else {
            cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
            cd b = -G[j1 + (long)i1 * ld + ld2]; if (i1 == j1) b = b + one;
            cd zxy = a * G[i1 + (long)j1 * ld + ld2] + b * G[i1 + (long)j1 * ld];
            cd zi = (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
            cd zj = (one - G[j1 + (long)j1 * ld + ld2]) - (one - G[j1 + (long)j1 * ld]);
            cd zz = a * G[i1 + (long)j1 * ld] + b * G[i1 + (long)j1 * ld + ld2] + zi * zj;
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

#endif  // ALFSB_WIP_CORR

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
    else
      oneapi::mkl::blas::column_major::gemm_batch(*Q, t(ta), t(tb), m, n, k, cd(ar, ai), A, lda, sA, B, ldb, sB, cd(br, bi), C, ldc,
                                                  sC, W);
  });
}
}
