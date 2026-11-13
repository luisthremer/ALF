// Intel batched study, Step 1 micro-benchmarks (one tile, one process; run several copies for sharing).
// Modes (argv[1]):
//   triad   : complex-double triad c = a + s*b over n elements, sizes from L2-resident (one G) to HBM-resident
//   gemm    : zgemm N x N x K (K = N: full; K = 16/32/64: rank-k flush), single and gemm_batch (strided) B = 1..16
//   launch  : empty SYCL kernel: submit+wait per kernel (sync) vs submit R kernels then one wait (async)
//   block   : emulation of one submatrix-update block (alf_la_sycl_blk.cpp sequence) for W walkers,
//             (a) today: W x the single-walker sequence, one D->H + wait per walker-block,
//             (b) batched: every op takes a walker dimension, one launch per op, one D->H + wait for all W
// Output: one line per measurement, "key value ..." (parsed by ib_sum.py).
#include <sycl/sycl.hpp>
#include <oneapi/mkl/blas.hpp>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
using cd = std::complex<double>;
using clk = std::chrono::steady_clock;
static double now() { return std::chrono::duration<double>(clk::now().time_since_epoch()).count(); }
namespace mb = oneapi::mkl::blas::column_major;
using oneapi::mkl::transpose;

static void fill(sycl::queue &q, cd *p, size_t n, double s) {
  q.parallel_for(sycl::range<1>(n), [=](sycl::id<1> i) {
     double x = (double)((i[0] * 2654435761ul) % 1000003ul) / 1000003.0;
     p[i] = cd(s * (x - 0.5), s * (0.25 - x * x));
   }).wait();
}

static void triad(sycl::queue &q) {
  size_t sizes[] = {576ul * 576, 2 * 576ul * 576, 1024ul * 1024, 8 * 576ul * 576, 16 * 1024ul * 1024, 64 * 1024ul * 1024};
  for (size_t n : sizes) {
    cd *a = sycl::malloc_device<cd>(n, q), *b = sycl::malloc_device<cd>(n, q), *c = sycl::malloc_device<cd>(n, q);
    fill(q, a, n, 1.0); fill(q, b, n, 2.0); fill(q, c, n, 0.0);
    cd s(0.5, -0.25);
    int reps = (int)std::max<size_t>(20, (size_t)4e9 / (n * 48));
    auto run = [&] { q.parallel_for(sycl::range<1>(n), [=](sycl::id<1> i) { c[i] = a[i] + s * b[i]; }); };
    for (int r = 0; r < 5; ++r) run();
    q.wait();
    double t0 = now();
    for (int r = 0; r < reps; ++r) run();
    q.wait();
    double t = (now() - t0) / reps;
    std::printf("triad n=%zu MB_per_array=%.1f t_us=%.2f GBs=%.1f\n", n, n * 16 / 1e6, t * 1e6, 3.0 * 16 * n / t / 1e9);
    sycl::free(a, q); sycl::free(b, q); sycl::free(c, q);
  }
}

static void gemm(sycl::queue &q, double tmin) {
  int Ns[] = {576, 1024};
  int Ks[] = {0, 16, 32, 64};
  int Bs[] = {1, 2, 4, 8, 16};
  for (int N : Ns)
    for (int K0 : Ks) {
      int K = K0 ? K0 : N;
      size_t sa = (size_t)N * K, sb = (size_t)K * N, sc = (size_t)N * N;
      int Bmax = 16;
      cd *A = sycl::malloc_device<cd>(sa * Bmax, q), *Bm = sycl::malloc_device<cd>(sb * Bmax, q), *C = sycl::malloc_device<cd>(sc * Bmax, q);
      fill(q, A, sa * Bmax, 1.0); fill(q, Bm, sb * Bmax, 1.0); fill(q, C, sc * Bmax, 1.0);
      cd al(-1.0, 0.0), be(1.0, 0.0);
      for (int B : Bs) {
        // (i) B separate gemm calls, async, one wait;  (ii) one gemm_batch strided call
        for (int mode = 0; mode < 2; ++mode) {
          if (B == 1 && mode == 1) continue;
          auto one = [&] {
            if (mode == 0)
              for (int b = 0; b < B; ++b)
                mb::gemm(q, transpose::nontrans, transpose::nontrans, N, N, K, al, A + sa * b, N, Bm + sb * b, K, be, C + sc * b, N);
            else
              mb::gemm_batch(q, transpose::nontrans, transpose::nontrans, N, N, K, al, A, N, sa, Bm, K, sb, be, C, N, sc, B);
          };
          for (int r = 0; r < 3; ++r) one();
          q.wait();
          int reps = 0;
          double t0 = now(), t;
          do { one(); ++reps; if (reps % 8 == 0) q.wait(); t = now() - t0; } while (t < tmin || reps % 8);
          q.wait();
          t = (now() - t0) / reps;
          double fl = 8.0 * N * N * (double)K * B;          // complex: 8 real flops per multiply-add
          double by = 16.0 * B * (sa + sb + 2 * sc);       // minimal traffic: read A,B,C, write C
          std::printf("gemm N=%d K=%d B=%d mode=%s t_us=%.1f us_per_prob=%.2f TFs=%.2f GBs=%.0f\n", N, K, B, mode ? "batch" : "loop", t * 1e6,
                      t * 1e6 / B, fl / t / 1e12, by / t / 1e9);
        }
      }
      sycl::free(A, q); sycl::free(Bm, q); sycl::free(C, q);
    }
}

static void launch(sycl::queue &q) {
  int R = 20000;
  int *d = sycl::malloc_device<int>(64, q);
  auto k = [&] { q.parallel_for(sycl::range<1>(64), [=](sycl::id<1> i) { d[i] = (int)i[0]; }); };
  for (int r = 0; r < 100; ++r) k();
  q.wait();
  double t0 = now();
  for (int r = 0; r < R; ++r) { k(); q.wait(); }
  double ts = (now() - t0) / R;
  t0 = now();
  for (int r = 0; r < R; ++r) k();
  double tsub = (now() - t0) / R;
  q.wait();
  double ta = (now() - t0) / R;
  // tiny H->D copy (pinned host) as the block path does, async and sync
  cd *h = sycl::malloc_host<cd>(64, q);
  cd *dd = sycl::malloc_device<cd>(64, q);
  t0 = now();
  for (int r = 0; r < R; ++r) { q.memcpy(dd, h, 64 * 16); q.wait(); }
  double tcs = (now() - t0) / R;
  t0 = now();
  for (int r = 0; r < R; ++r) q.memcpy(dd, h, 64 * 16);
  q.wait();
  double tca = (now() - t0) / R;
  std::printf("launch sync_us=%.2f async_submit_us=%.2f async_total_us=%.2f copy_sync_us=%.2f copy_async_us=%.2f\n", ts * 1e6, tsub * 1e6,
              ta * 1e6, tcs * 1e6, tca * 1e6);
  sycl::free(d, q); sycl::free(h, q); sycl::free(dd, q);
}

// ---- update-block emulation -------------------------------------------------------------------------------------
// Per walker: G (N x N x 2 flavours), block of b sites, flush of k columns per flavour (here: every block flushes k).
// Sequence today per walker-block (alf_la_sycl_blk.cpp, p7 NOTES section 4): 4 small H->D uploads, gather (2 x b x b),
// D->H of the gathered block + wait (host Metropolis), scale (2 x b x N), fgather x2 (k x N), 4 rank-k gemms (2 per flavour,
// m = n = N/2.. here N x N x k/2 to keep flops of 2 gemms per flavour), copies.
struct Blk {
  int N, b, k, W;
  cd *G, *gbk, *cbk, *Cc, *Qt, *hbuf, *hup, *Sx, *Mk;
  int *bsk, *S, *slot;
};
static void blk_alloc(sycl::queue &q, Blk &x) {
  size_t NN = (size_t)x.N * x.N;
  x.G = sycl::malloc_device<cd>(NN * 2 * x.W, q);
  x.gbk = sycl::malloc_device<cd>((size_t)2 * x.b * x.b * x.W, q);
  x.cbk = sycl::malloc_device<cd>((size_t)2 * 2 * x.b * x.W, q);
  x.Cc = sycl::malloc_device<cd>((size_t)x.N * x.k * 2 * x.W, q);
  x.Qt = sycl::malloc_device<cd>((size_t)x.N * x.k * 2 * x.W, q);
  x.Sx = sycl::malloc_device<cd>((size_t)x.N * x.k * 2 * x.W, q);
  x.Mk = sycl::malloc_device<cd>((size_t)x.k * x.k * 2 * x.W, q);
  fill(q, x.Mk, (size_t)x.k * x.k * 2 * x.W, 0.01);
  x.bsk = sycl::malloc_device<int>((size_t)2 * x.b * x.W, q);
  x.S = sycl::malloc_device<int>((size_t)2 * x.k * x.W, q);
  x.slot = sycl::malloc_device<int>((size_t)x.N * 2 * x.W, q);
  x.hbuf = sycl::malloc_host<cd>((size_t)2 * x.b * x.b * x.W, q);
  x.hup = sycl::malloc_host<cd>((size_t)2 * 2 * x.b * x.W + 4096, q);
  fill(q, x.G, NN * 2 * x.W, 0.01);
  fill(q, x.Cc, (size_t)x.N * x.k * 2 * x.W, 0.01);
  fill(q, x.Qt, (size_t)x.N * x.k * 2 * x.W, 0.01);
  std::vector<int> hb(2 * x.b * x.W), hs(2 * x.k * x.W);
  for (size_t i = 0; i < hb.size(); ++i) hb[i] = 1 + (int)((i * 37) % x.N);
  for (size_t i = 0; i < hs.size(); ++i) hs[i] = 1 + (int)((i * 53) % x.N);
  q.memcpy(x.bsk, hb.data(), hb.size() * 4).wait();
  q.memcpy(x.S, hs.data(), hs.size() * 4).wait();
  for (int i = 0; i < 2 * 2 * x.b * x.W; ++i) x.hup[i] = cd(1.0 + 1e-3 * i, 0.0);
  q.memcpy(x.cbk, x.hup, (size_t)2 * 2 * x.b * x.W * 16).wait();
}
// one op = all walkers w in [w0, w0+nw) in one launch
static void blk_ops(sycl::queue &q, Blk &x, int w0, int nw, bool wait_dh) {
  const int N = x.N, b = x.b, k = x.k;
  const size_t NN = (size_t)N * N;
  cd *G = x.G, *gbk = x.gbk, *cbk = x.cbk, *Cc = x.Cc, *Qt = x.Qt;
  int *bsk = x.bsk, *S = x.S, *slot = x.slot;
  // 4 uploads (sizes as today: small)
  for (int u = 0; u < 4; ++u) q.memcpy(cbk + (size_t)2 * 2 * b * w0, x.hup + (size_t)2 * 2 * b * w0, (size_t)2 * 2 * b * nw * 16 / 4);
  // gather
  q.parallel_for(sycl::range<3>(nw * 2, b, b), [=](sycl::id<3> id) {
    long wf = id[0] + 2l * w0, c = id[1], a = id[2];
    gbk[a + (long)b * (c + (long)b * wf)] = G[(bsk[a + (long)b * wf] - 1) + (long)N * (bsk[c + (long)b * wf] - 1) + NN * wf];
  });
  q.memcpy(x.hbuf + (size_t)2 * b * b * w0, gbk + (size_t)2 * b * b * w0, (size_t)2 * b * b * nw * 16);
  if (wait_dh) q.wait();  // host Metropolis needs the block
  // scale
  q.parallel_for(sycl::range<3>(nw * 2, b, N), [=](sycl::id<3> id) {
    long wf = id[0] + 2l * w0, a = id[1], j = id[2];
    long sa = bsk[a + (long)b * wf] - 1;
    cd *Gf = G + NN * wf;
    Gf[sa + (long)N * j] = cbk[a + 2l * b * wf] * Gf[sa + (long)N * j];
  });
  // flush gathers (2 per walker as today: one per flavour)
  for (int f = 0; f < 2; ++f)
    q.parallel_for(sycl::range<2>(nw * k, N), [=](sycl::id<2> id) {
      long w = id[0] / k + w0, m = id[0] % k, i = id[1];
      long wf = 2 * w + f;
      long s = S[m + (long)k * wf];
      const cd *Gf = G + NN * wf;
      Cc[i + (long)N * (m + (long)k * wf)] = Gf[i + (long)N * (s - 1)];
      cd qv = -Gf[(s - 1) + (long)N * i];
      if (i + 1 == s) { qv += 1.0; slot[i + (long)N * wf] = 0; }
      Qt[i + (long)N * (m + (long)k * wf)] = qv;
    });
  cd al(-1e-3, 0.0), be(1.0, 0.0), one(1.0, 0.0), zero(0.0, 0.0);
  if (getenv("IB_BLK_REAL")) {
    // as gpu_lowrank_mod ab_flush_dev: per flavour gemm1 Sx(N x k) = Qt(N x k) * M(k x k)^T, gemm2 G -= Cc * Sx^T (one G pass)
    cd *Sx = x.Sx, *Mk = x.Mk;
    if (nw == 1) {
      for (int f = 0; f < 2; ++f) {
        long wf = 2l * w0 + f;
        mb::gemm(q, transpose::nontrans, transpose::trans, N, k, k, one, Qt + (size_t)N * k * wf, N, Mk + (size_t)k * k * wf, k, zero, Sx + (size_t)N * k * wf, N);
        mb::gemm(q, transpose::nontrans, transpose::trans, N, N, k, al, Cc + (size_t)N * k * wf, N, Sx + (size_t)N * k * wf, N, be, G + NN * wf, N);
      }
    } else {
      mb::gemm_batch(q, transpose::nontrans, transpose::trans, N, k, k, one, Qt + (size_t)N * k * 2 * w0, N, (int64_t)N * k, Mk + (size_t)k * k * 2 * w0, k,
                     (int64_t)k * k, zero, Sx + (size_t)N * k * 2 * w0, N, (int64_t)N * k, 2 * nw);
      mb::gemm_batch(q, transpose::nontrans, transpose::trans, N, N, k, al, Cc + (size_t)N * k * 2 * w0, N, (int64_t)N * k, Sx + (size_t)N * k * 2 * w0, N,
                     (int64_t)N * k, be, G + NN * 2 * w0, N, (int64_t)NN, 2 * nw);
    }
  } else {
  // variant A (job 560867): 2 gemms per flavour of rank k/2, each a full G pass (4 G passes per walker-block: pessimistic)
  int kh = k / 2;
  if (nw == 1) {
    for (int f = 0; f < 2; ++f)
      for (int h = 0; h < 2; ++h) {
        long wf = 2l * w0 + f;
        mb::gemm(q, transpose::nontrans, transpose::trans, N, N, kh, al, Cc + (size_t)N * (k * wf + h * kh), N, Qt + (size_t)N * (k * wf + h * kh), N, be,
                 G + NN * wf, N);
      }
  } else {
    for (int h = 0; h < 2; ++h)
      mb::gemm_batch(q, transpose::nontrans, transpose::trans, N, N, kh, al, Cc + (size_t)N * (k * 2l * w0 + h * kh), N, (int64_t)N * k,
                     Qt + (size_t)N * (k * 2l * w0 + h * kh), N, (int64_t)N * k, be, G + NN * 2 * w0, N, (int64_t)NN, 2 * nw);
  }
  }
  // trailing copies (today: several small D->H/H->D copies per block)
  for (int u = 0; u < 4; ++u) q.memcpy(x.hup + 4096 - 64, cbk + (size_t)2 * 2 * b * w0, 64 * 16 / 4);
}
static void block(sycl::queue &q, int N, int b, int k, double tmin) {
  int Ws[] = {1, 2, 4, 8, 16};
  for (int W : Ws) {
    Blk x{N, b, k, W};
    blk_alloc(q, x);
    for (int mode = 0; mode < 2; ++mode) {
      auto one = [&] {
        if (mode == 0)
          for (int w = 0; w < W; ++w) blk_ops(q, x, w, 1, true);
        else
          blk_ops(q, x, 0, W, true);
      };
      for (int r = 0; r < 5; ++r) one();
      q.wait();
      int reps = 0;
      double t0 = now(), t;
      do { one(); ++reps; t = now() - t0; } while (t < tmin);
      q.wait();
      t = (now() - t0) / reps;
      std::printf("block N=%d b=%d k=%d W=%d mode=%s t_us=%.1f us_per_walker_block=%.2f\n", N, b, k, W, mode ? "batched" : "today", t * 1e6, t * 1e6 / W);
    }
    // pure device time of the batched op set (no host waits inside): submit many, wait once
    {
      int reps = 0;
      double t0 = now(), t;
      do { blk_ops(q, x, 0, W, false); ++reps; t = now() - t0; } while (t < tmin);
      q.wait();
      t = (now() - t0) / reps;
      std::printf("block N=%d b=%d k=%d W=%d mode=batched_nowait t_us=%.1f us_per_walker_block=%.2f\n", N, b, k, W, t * 1e6, t * 1e6 / W);
    }
    sycl::free(x.G, q); sycl::free(x.gbk, q); sycl::free(x.cbk, q); sycl::free(x.Cc, q); sycl::free(x.Qt, q); sycl::free(x.Sx, q); sycl::free(x.Mk, q);
    sycl::free(x.bsk, q); sycl::free(x.S, q); sycl::free(x.slot, q); sycl::free(x.hbuf, q); sycl::free(x.hup, q);
  }
}


// ---- observables corr kernel (SpinMz kind of gpu_obs_mod corr): today's mapping vs a column-coalesced two-pass form ----------
// square lattice L x L, N = L^2, Norb = 1. acc(d) += sum_I f(G(J,I,:), G(I,J,:), G(I,I,:), G(J,J,:)), J = I + d.
// A (today): one work-group per output d, work-items over I (both G(J,I) and G(I,J) uncoalesced).
// B: B1 one work-item per (J fastest, I): reads G(J,I) contiguous, writes T(d, I) at [I*N + d] (d = J - I: contiguous in J except
//    at row wraps); B2 reduces over I in C chunks with work-items over d (coalesced), B3 adds the C partials.
static void corrbench(sycl::queue &q, int L, double tmin) {
  const int N = L * L;
  const size_t NN = (size_t)N * N;
  cd *G = sycl::malloc_device<cd>(NN * 2, q);
  fill(q, G, NN * 2, 0.01);
  std::vector<int> hjc(NN), hinv(NN);
  for (int I = 0; I < N; ++I)
    for (int d = 0; d < N; ++d) {
      int x = I % L, y = I / L, dx = d % L, dy = d / L;
      int J = (x + dx) % L + L * ((y + dy) % L);
      hjc[I + (size_t)N * d] = J; hinv[J + (size_t)N * I] = d;
    }
  int *jc = sycl::malloc_device<int>(NN, q), *inv = sycl::malloc_device<int>(NN, q);
  q.memcpy(jc, hjc.data(), NN * 4); q.memcpy(inv, hinv.data(), NN * 4); q.wait();
  cd *acc = sycl::malloc_device<cd>((size_t)3 * N, q);
  const int C = 16;
  cd *T = sycl::malloc_device<cd>(NN * 3, q), *P = sycl::malloc_device<cd>((size_t)3 * N * C, q);
  auto term = [](const cd *G, size_t NN, int N, int I, int J, cd &s1, cd &s2, cd &s3) {
    cd a = -G[J + (size_t)N * I], b = -G[J + (size_t)N * I + NN];
    if (I == J) { a += 1.0; b += 1.0; }
    cd gij1 = G[I + (size_t)N * J], gij2 = G[I + (size_t)N * J + NN];
    cd zxy = a * gij2 + b * gij1;
    cd zi = (1.0 - G[I + (size_t)N * I + NN]) - (1.0 - G[I + (size_t)N * I]);
    cd zj = (1.0 - G[J + (size_t)N * J + NN]) - (1.0 - G[J + (size_t)N * J]);
    cd zz = a * gij1 + b * gij2 + zi * zj;
    s1 = zz; s2 = zxy; s3 = (2.0 * zxy + zz) / 3.0;
  };
  const int WG = 256;
  auto runA = [&] {
    q.parallel_for(sycl::nd_range<1>((size_t)N * WG, WG), [=](sycl::nd_item<1> it) {
      int d = it.get_group(0);
      double r[6] = {0, 0, 0, 0, 0, 0};
      for (int I = it.get_local_id(0); I < N; I += WG) {
        int J = jc[I + (size_t)N * d];
        cd s1, s2, s3;
        term(G, NN, N, I, J, s1, s2, s3);
        r[0] += s1.real(); r[1] += s1.imag(); r[2] += s2.real(); r[3] += s2.imag(); r[4] += s3.real(); r[5] += s3.imag();
      }
      for (int c = 0; c < 6; ++c) r[c] = sycl::reduce_over_group(it.get_group(), r[c], sycl::plus<double>());
      if (it.get_local_id(0) == 0) {
        acc[d] += cd(r[0], r[1]); acc[d + N] += cd(r[2], r[3]); acc[d + 2 * N] += cd(r[4], r[5]);
      }
    });
  };
  auto runB = [&] {
    q.parallel_for(sycl::range<2>(N, N), [=](sycl::id<2> id) {
      int I = id[0], J = id[1];
      int d = inv[J + (size_t)N * I];
      cd s1, s2, s3;
      term(G, NN, N, I, J, s1, s2, s3);
      size_t o = (size_t)I * N + d;
      T[o] = s1; T[o + NN] = s2; T[o + 2 * NN] = s3;
    });
    q.parallel_for(sycl::range<2>(C, N), [=](sycl::id<2> id) {
      int c = id[0], d = id[1];
      int i0 = (int)((long)N * c / C), i1 = (int)((long)N * (c + 1) / C);
      cd r1 = 0, r2 = 0, r3 = 0;
      for (int I = i0; I < i1; ++I) { size_t o = (size_t)I * N + d; r1 += T[o]; r2 += T[o + NN]; r3 += T[o + 2 * NN]; }
      P[d + (size_t)N * (c + C * 0)] = r1; P[d + (size_t)N * (c + C * 1)] = r2; P[d + (size_t)N * (c + C * 2)] = r3;
    });
    q.parallel_for(sycl::range<2>(3, N), [=](sycl::id<2> id) {
      int k = id[0], d = id[1];
      cd r = 0;
      for (int c = 0; c < C; ++c) r += P[d + (size_t)N * (c + C * k)];
      acc[d + (size_t)N * k] += r;
    });
  };
  std::vector<cd> ra(3 * N), rb(3 * N);
  for (int v = 0; v < 2; ++v) {
    q.memset(acc, 0, 16 * 3 * N).wait();
    if (v == 0) runA(); else runB();
    q.wait();
    q.memcpy(v == 0 ? ra.data() : rb.data(), acc, 16 * 3 * N).wait();
    for (int r = 0; r < 3; ++r) { if (v == 0) runA(); else runB(); }
    q.wait();
    int reps = 0;
    double t0 = now(), t;
    do { if (v == 0) runA(); else runB(); ++reps; if (reps % 16 == 0) q.wait(); t = now() - t0; } while (t < tmin || reps % 16);
    q.wait();
    t = (now() - t0) / reps;
    std::printf("corr L=%d N=%d variant=%s t_us=%.1f\n", L, N, v == 0 ? "A_today" : "B_twopass", t * 1e6);
  }
  double md = 0, mx = 0;
  for (int i = 0; i < 3 * N; ++i) { md = std::max(md, std::abs(ra[i] - rb[i])); mx = std::max(mx, std::abs(ra[i])); }
  std::printf("corr L=%d check maxabs_diff=%.3e max=%.3e\n", L, md, mx);
  sycl::free(G, q); sycl::free(jc, q); sycl::free(inv, q); sycl::free(acc, q); sycl::free(T, q); sycl::free(P, q);
}

int main(int argc, char **argv) {
  std::string mode = argc > 1 ? argv[1] : "all";
  double tmin = argc > 2 ? atof(argv[2]) : 1.0;
#ifdef IB_OMPINIT
  // oneAPI 2025.3 with the UR Level-Zero v1 adapter: a pure SYCL program fails at platform init (UR_RESULT_ERROR_UNINITIALIZED);
  // ALF initialises the device through libomptarget first. Do the same.
  {
    int x = 0;
#pragma omp target map(tofrom : x)
    x = 1;
  }
#endif
  sycl::queue q(sycl::gpu_selector_v, sycl::property::queue::in_order());
  auto dev = q.get_device();
  std::printf("device %s CUs=%u\n", dev.get_info<sycl::info::device::name>().c_str(), dev.get_info<sycl::info::device::max_compute_units>());
  std::fflush(stdout);
  if (mode == "triad" || mode == "all") triad(q);
  std::fflush(stdout);
  if (mode == "launch" || mode == "all") launch(q);
  std::fflush(stdout);
  if (mode == "gemm" || mode == "all") gemm(q, tmin);
  std::fflush(stdout);
  if (mode == "corr") {
    corrbench(q, 24, tmin);
    corrbench(q, 32, tmin);
  }
  if (mode == "block" || mode == "all") {
    int kb = getenv("IB_BLK_K") ? atoi(getenv("IB_BLK_K")) : 32;
    block(q, 576, 32, kb, tmin);
    block(q, 1024, 32, kb, tmin);
  }
  return 0;
}
