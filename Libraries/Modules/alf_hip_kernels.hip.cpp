// ALF GPU=amd: small device kernels written in HIP instead of OpenMP target regions.
// A libgomp-GCN target region with a parallel loop costs ~150 us to launch; a HIP kernel launch + hipDeviceSynchronize ~15 us.
// Every function takes device pointers (from use_device_addr in alf_hipk_mod) and finishes (synchronises) before it returns,
// because the OpenMP target regions run on their own queues. Index conventions follow the Fortran arrays (column major).
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>

struct cd {
  double r, i;
  __device__ cd() : r(0), i(0) {}
  __device__ cd(double a, double b) : r(a), i(b) {}
};
__device__ inline cd operator+(cd a, cd b) { return cd(a.r + b.r, a.i + b.i); }
__device__ inline cd operator-(cd a, cd b) { return cd(a.r - b.r, a.i - b.i); }
__device__ inline cd operator-(cd a) { return cd(-a.r, -a.i); }
__device__ inline cd operator*(cd a, cd b) { return cd(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
__device__ inline cd operator*(double s, cd a) { return cd(s * a.r, s * a.i); }
__device__ inline cd operator*(cd a, double s) { return cd(s * a.r, s * a.i); }
__device__ inline cd operator/(cd a, double s) { return cd(a.r / s, a.i / s); }
__device__ inline cd& operator+=(cd& a, cd b) { a.r += b.r; a.i += b.i; return a; }

// ---- deferred synchronisation -------------------------------------------------------------------------------------
// All HIP kernels, hipBLAS calls and the QR run on stream 0 and are therefore ordered among themselves. Only the OpenMP
// target regions (libgomp queue) and the copies of the libgomp plugin need the HIP work to be finished: the patched GCN
// plugin calls alf_hip_presync() before every OpenMP kernel launch and before every device copy/free
// (alf_gcn_set_presync, registered in alfhk_register_presync). So the wrappers only mark the work as pending.
// ALF_HIP_ASYNC=0 restores a synchronisation after every call (debugging, A/B timing).
extern "C" {
int alf_hip_pending = 0;
static int g_async = -1;
static int g_registered = 0;
int alf_hip_presync() {
  if (!alf_hip_pending) return 0;
  alf_hip_pending = 0;
  hipError_t e = hipDeviceSynchronize();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: sync error %s\n", hipGetErrorString(e)); return 2; }
  return 0;
}
int alf_hip_fin() {
  hipError_t e = hipGetLastError();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: launch error %s\n", hipGetErrorString(e)); return 1; }
  if (g_async < 0) { const char* v = getenv("ALF_HIP_ASYNC"); g_async = (v && v[0] == '0') ? 0 : g_registered;   // without the plugin hook every call must synchronise
    printf(" alf_la: HIP deferred synchronisation: %s\n", g_async ? "ON" : (g_registered ? "OFF (ALF_HIP_ASYNC=0)" : "OFF (libgomp plugin has no presync hook)")); }
  if (g_async) { alf_hip_pending = 1; return 0; }
  e = hipDeviceSynchronize();
  if (e != hipSuccess) { fprintf(stderr, "alf_hip_kernels: sync error %s\n", hipGetErrorString(e)); return 2; }
  return 0;
}
// hook the presync into the (patched) libgomp GCN plugin; returns 1 if the plugin has the hook, else 0 (then every call syncs)
int alfhk_register_presync() {
  void* h = dlopen("libgomp-plugin-gcn.so.1", RTLD_NOLOAD | RTLD_LAZY);
  if (!h) return 0;
  void (*set)(int (*)()) = (void (*)(int (*)()))dlsym(h, "alf_gcn_set_presync");
  if (!set) return 0;
  set(alf_hip_presync);
  g_registered = 1;
  return 1;
}
}
static int fin() { return alf_hip_fin(); }
static inline unsigned nb(long n, int t) { return (unsigned)((n + t - 1) / t); }
#define LAUNCH(k, grid, blk, ...) hipLaunchKernelGGL(k, dim3(grid), dim3(blk), 0, 0, __VA_ARGS__)

// ---------------------------------------------------------------------------------------------------------------------
__global__ void k_copyz(long n, const cd* a, cd* b) {
  long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) b[i] = a[i];
}

// gbd(a,c,fl) = G(bsd(a,fl), bsd(c,fl), nf(fl)),  gbd(nd,nd,2), bsd(nd,2), G(ld,ld,*)
__global__ void k_blk_gather(int nfl, int nfa, int nfb, int b, int nd, int ld, const cd* G, const int* bsd, cd* gbd) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  long tot = (long)nfl * b * b;
  if (t >= tot) return;
  int a = t % b, c = (t / b) % b, fl = t / ((long)b * b);
  int nf = (fl == 1) ? nfb : nfa;
  int ra = bsd[a + fl * nd] - 1, rc = bsd[c + fl * nd] - 1;
  gbd[a + (long)c * nd + (long)fl * nd * nd] = G[ra + (long)rc * ld + (long)(nf - 1) * ld * ld];
}

// G(bsd(a,fl), j, nf) *= cbuf(a,fl)   (rows),  G(i, bsd(a,fl), nf) *= cbuf(nd+a,fl)  (columns);  cbuf(cbld, 2)
__global__ void k_blk_scale_rows(int nf, int fl, int ld, int b, int nd, int cbld, cd* G, const int* bsd, const cd* cbuf) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)b * ld) return;
  int a = t % b, j = t / b;
  int r = bsd[a + (fl - 1) * nd] - 1;
  cd f = cbuf[a + (long)(fl - 1) * cbld];
  cd* p = G + r + (long)j * ld + (long)(nf - 1) * ld * ld;
  *p = f * (*p);
}
__global__ void k_blk_scale_cols(int nf, int fl, int ld, int b, int nd, int cbld, cd* G, const int* bsd, const cd* cbuf) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)b * ld) return;
  int i = t % ld, a = t / ld;
  int c = bsd[a + (fl - 1) * nd] - 1;
  cd f = cbuf[nd + a + (long)(fl - 1) * cbld];
  cd* p = G + i + (long)c * ld + (long)(nf - 1) * ld * ld;
  *p = f * (*p);
}

// Fused block end for one flavour: rows of the block sites scaled by cbuf(a), then columns by cbuf(nd+a) (same operation order as
// k_blk_scale_rows followed by k_blk_scale_cols, bit-identical), and the flush gathers C(i,m) = G(i,S(m)), Qt(i,m) = -G(S(m),i) (+1,
// slot = 0 on the diagonal) from the final values.  The pending sites S are a subset of the block sites, so every element needed by the
// gather is written here.  Threads [0,b*ld): rows part; [b*ld,2*b*ld): columns part for rows outside the block.
__global__ void k_blk_end(int nf, int fl, int ld, int b, int nd, int cbld, int kk, int nslot_ld, cd* G, const int* bsd, const cd* cbuf,
                          const int* S, cd* C, cd* Qt, int* slot) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  long tot = (long)b * ld;
  const int* bs = bsd + (long)(fl - 1) * nd;
  const cd* cb = cbuf + (long)(fl - 1) * cbld;
  cd* Gf = G + (long)(nf - 1) * ld * ld;
  if (t < tot) {
    int a = t % b, j = t / b;
    int r = bs[a] - 1;
    cd v = cb[a] * Gf[r + (long)j * ld];
    for (int q = 0; q < b; ++q) if (bs[q] - 1 == j) { v = cb[nd + q] * v; break; }
    Gf[r + (long)j * ld] = v;
    for (int m = 0; m < kk; ++m) {
      int s = S[m] - 1;
      if (s == r) { cd w = -v; if (j == r) { w.r += 1.0; slot[j + (long)(nf - 1) * nslot_ld] = 0; } Qt[j + (long)m * ld] = w; }
      if (s == j) C[r + (long)m * ld] = v;
    }
  } else if (t < 2 * tot) {
    long t2 = t - tot;
    int i = t2 % ld, a = t2 / ld;
    for (int q = 0; q < b; ++q) if (bs[q] - 1 == i) return;
    int c = bs[a] - 1;
    cd v = cb[nd + a] * Gf[i + (long)c * ld];
    Gf[i + (long)c * ld] = v;
    for (int m = 0; m < kk; ++m) if (S[m] - 1 == c) C[i + (long)m * ld] = v;
  }
}

// C(i,m) = G(i,S(m)); Qt(i,m) = -G(S(m),i) (+1 if i==S(m), and slot(i,nf) = 0);  C, Qt (ld,*), G (ld,*), slot (nslot_ld, *)
__global__ void k_flush_gather(int ld, int nrow, int k, int nf, int nslot_ld, const cd* G, const int* S, cd* C, cd* Qt, int* slot) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)k * nrow) return;
  int i = t % nrow, m = t / nrow;
  int s = S[m] - 1;
  C[i + (long)m * ld] = G[i + (long)s * ld];
  cd q = -G[s + (long)i * ld];
  if (i == s) { q.r += 1.0; slot[i + (long)(nf - 1) * nslot_ld] = 0; }
  Qt[i + (long)m * ld] = q;
}

// M(P(l), j) = W(l, j)
__global__ void k_scatter_rows(int nd, int ncol, const int* P, cd* M, int ld, const cd* W, int ldw) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)nd * ncol) return;
  int l = t % nd, j = t / nd;
  M[(P[l] - 1) + (long)j * ld] = W[l + (long)j * ldw];
}

// ---- observables --------------------------------------------------------------------------------------------------
#define RB 128
// One block per output element m = (imj,no,no2) (gather over the lattice sites I) for m < MT; blocks MT..MT+Norb-1 do the acc0 sums (kinds 3, 4).
__global__ void k_corr(int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld,
                       const cd* G, const int* jc, const int* inv, cd* acc, cd* acc0) {
  __shared__ double red[6][RB];
  long MT = (long)N * Norb * Norb;
  long m = blockIdx.x;
  bool zero_part = (m >= MT);
  int imj, no, no2;
  if (!zero_part) { imj = m % N; no = (m / N) % Norb; no2 = m / ((long)N * Norb); }
  else { no = (int)(m - MT); imj = 0; no2 = 0; }
  cd s1, s2, s3;
  const long ld2 = (long)ld * ld;
  for (int I = threadIdx.x; I < N; I += blockDim.x) {
    if (!zero_part) {
      int J = jc[I + (long)imj * N] - 1;
      int i1 = inv[I + (long)no * N] - 1, j1 = inv[J + (long)no2 * N] - 1;
      cd one(1.0, 0.0);
      if (kind == 1) {
        cd z;
        for (int nf = 0; nf < NFL; nf++) { cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one; z += a; }
        s1 += z * ns;
      } else if (kind == 2) {
        cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
        s1 += (a * G[i1 + (long)j1 * ld]) * ns;
      } else if (kind == 3) {
        cd zi, zj, z;
        for (int nf = 0; nf < NFL; nf++) {
          zi += one - G[i1 + (long)i1 * ld + nf * ld2];
          zj += one - G[j1 + (long)j1 * ld + nf * ld2];
          cd a = -G[j1 + (long)i1 * ld + nf * ld2]; if (i1 == j1) a = a + one;
          z += a * G[i1 + (long)j1 * ld + nf * ld2];
        }
        s1 += ((zi * ns) * (zj * ns)) + z * ns;
      } else {
        cd a = -G[j1 + (long)i1 * ld]; if (i1 == j1) a = a + one;
        cd b = -G[j1 + (long)i1 * ld + ld2]; if (i1 == j1) b = b + one;
        cd zxy = a * G[i1 + (long)j1 * ld + ld2] + b * G[i1 + (long)j1 * ld];
        cd zi = (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
        cd zj = (one - G[j1 + (long)j1 * ld + ld2]) - (one - G[j1 + (long)j1 * ld]);
        cd zz = a * G[i1 + (long)j1 * ld] + b * G[i1 + (long)j1 * ld + ld2] + zi * zj;
        s1 += zz; s2 += zxy; s3 += (2.0 * zxy + zz) / 3.0;
      }
    } else {
      int i1 = inv[I + (long)no * N] - 1;
      cd one(1.0, 0.0);
      if (kind == 3) {
        cd zi;
        for (int nf = 0; nf < NFL; nf++) zi += one - G[i1 + (long)i1 * ld + nf * ld2];
        s1 += zi * ns;
      } else {
        s1 += (one - G[i1 + (long)i1 * ld + ld2]) - (one - G[i1 + (long)i1 * ld]);
      }
    }
  }
  red[0][threadIdx.x] = s1.r; red[1][threadIdx.x] = s1.i;
  red[2][threadIdx.x] = s2.r; red[3][threadIdx.x] = s2.i;
  red[4][threadIdx.x] = s3.r; red[5][threadIdx.x] = s3.i;
  __syncthreads();
  for (int st = RB / 2; st > 0; st >>= 1) {
    if (threadIdx.x < st) for (int q = 0; q < 6; q++) red[q][threadIdx.x] += red[q][threadIdx.x + st];
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    cd w(wr, wi);
    cd r1(red[0][0], red[1][0]), r2(red[2][0], red[3][0]), r3(red[4][0], red[5][0]);
    if (!zero_part) {
      acc[m + (long)(k - 1) * MT] += r1 * w;
      if (kind == 4) { acc[m + (long)k * MT] += r2 * w; acc[m + (long)(k + 1) * MT] += r3 * w; }
    } else {
      acc0[no + (long)(k - 1) * Norb] += r1 * w;
      if (kind == 4) acc0[no + (long)(k + 1) * Norb] += (r1 * w) / sq3;
    }
  }
}

// Zkin (nk terms), Zpot, Zrho in one block; out[0..5] = re/im of Zkin(before N_SUN), Zpot, Zrho
__global__ void k_scalars(int nk, const int* hi, const int* hj, const int* hf, const cd* hc, int N, int Norb, int NFL, int dec,
                          int ld, const cd* G, const int* inv, const double* uo, double* out) {
  __shared__ double red[6][RB];
  const long ld2 = (long)ld * ld;
  cd zk, zp, zr;
  cd one(1.0, 0.0);
  for (int k = threadIdx.x; k < nk; k += blockDim.x) {
    int i1 = hi[k] - 1, j1 = hj[k] - 1, nf = hf[k] - 1;
    cd zc = -G[j1 + (long)i1 * ld + nf * ld2];
    if (i1 == j1) zc = zc + one;
    zk += hc[k] * zc;
  }
  for (int t = threadIdx.x; t < N * Norb; t += blockDim.x) {
    int I = t % N, no = t / N;
    int i1 = inv[I + (long)no * N] - 1;
    cd zc = one - G[i1 + (long)i1 * ld];
    cd zt = one - G[i1 + (long)i1 * ld + (long)(dec - 1) * ld2];
    zp += (zc * zt) * uo[no];
    for (int nf = 0; nf < NFL; nf++) zr += one - G[i1 + (long)i1 * ld + nf * ld2];
  }
  red[0][threadIdx.x] = zk.r; red[1][threadIdx.x] = zk.i;
  red[2][threadIdx.x] = zp.r; red[3][threadIdx.x] = zp.i;
  red[4][threadIdx.x] = zr.r; red[5][threadIdx.x] = zr.i;
  __syncthreads();
  for (int st = RB / 2; st > 0; st >>= 1) {
    if (threadIdx.x < st) for (int q = 0; q < 6; q++) red[q][threadIdx.x] += red[q][threadIdx.x + st];
    __syncthreads();
  }
  if (threadIdx.x == 0) for (int q = 0; q < 6; q++) out[q] = red[q][0];
}

// local observables: d = (G(i1,i1,dec) - G(i1,i1,1)) * w; lacc1(I,no) += d; lacc2(I,it,no) += wt*d
__global__ void k_local(int N, int Norb, int lNt, int it, int dec, double wt, double wr, double wi, int ld, const int* inv,
                        const cd* G, cd* lacc1, cd* lacc2) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)N * Norb) return;
  int I = t % N, no = t / N;
  int i1 = inv[I + (long)no * N] - 1;
  const long ld2 = (long)ld * ld;
  cd d = (G[i1 + (long)i1 * ld + (long)(dec - 1) * ld2] - G[i1 + (long)i1 * ld]) * cd(wr, wi);
  lacc1[I + (long)no * N] += d;
  lacc2[I + (long)(it - 1) * N + (long)no * N * lNt] += wt * d;
}


// ---- stabilisation helpers ----------------------------------------------------------------------------------------
__device__ inline cd conjg(cd a) { return cd(a.r, -a.i); }
__device__ inline double cabsd(cd a) { return hypot(a.r, a.i); }
__device__ inline cd cdiv(cd a, cd b) {   // Smith's algorithm
  if (fabs(b.r) >= fabs(b.i)) { double q = b.i / b.r, d = b.r + b.i * q; return cd((a.r + a.i * q) / d, (a.i - a.r * q) / d); }
  double q = b.r / b.i, d = b.r * q + b.i;
  return cd((a.r * q + a.i) / d, (a.i * q - a.r) / d);
}

// U(i,j) /= real(D(i)) for j >= i
__global__ void k_divrows(int n, int ld, cd* U, const cd* D) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  if (j >= i) U[i + (long)j * ld] = U[i + (long)j * ld] / D[i].r;
}
__global__ void k_rowdiv(int n, int lda, cd* A, const cd* D) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  A[i + (long)j * lda] = cdiv(A[i + (long)j * lda], D[i]);
}
__global__ void k_scale_cols(int n, int ld, cd* U, const cd* D) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  U[i + (long)j * ld] = U[i + (long)j * ld] * D[j];
}
__global__ void k_tpup(int n, cd* T, const cd* R, const cd* Dr, const cd* Dl) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  cd v = (Dr[i] * T[i + (long)j * n]) * Dl[j];
  T[i + (long)j * n] = v + R[i + (long)j * n];
}
__global__ void k_ns(int n, cd* A) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  cd v = -0.5 * A[i + (long)j * n];
  if (i == j) v.r += 1.5;
  A[i + (long)j * n] = v;
}
// complex product over the 64 lanes of a wave (butterfly; every lane returns a product, lane 0 is used)
__device__ inline cd wave_prod(cd x) {
  for (int o = 32; o > 0; o >>= 1) { cd y(__shfl_xor(x.r, o, 64), __shfl_xor(x.i, o, 64)); x = x * y; }
  return x;
}
// phase accumulation: cs(3) *= ...  The factors of the N reflectors are independent: one wave forms strided partial
// products and combines them (the order of the rounding differs from the serial CPU loop at the 1e-16 level)
__global__ void k_cgr_phase(int n, int ld, const cd* T, const cd* tau, cd* cs, int nvar, int sgn, int bad) {
  cd acc(1.0, 0.0);
  for (int i = threadIdx.x; i < n; i += 64) {
    cd z = tau[i];
    cd tii = T[i + (long)i * ld];
    double at = cabsd(tii);
    cd f;
    if (nvar == 1) f = tii / at;
    else { f = conjg(tii) / at; z = conjg(z); }
    if (z.r != 0.0 || z.i != 0.0) {
      double x = cabsd(z);
      cd q = z / x;
      cd w = cd(1.0, 0.0) - 2.0 * (q * (z.r / x));
      f = f * (w / cabsd(w));
    }
    acc = acc * f;
  }
  acc = wave_prod(acc);
  if (threadIdx.x == 0) {
    cd ph = cs[2];
    if (sgn < 0) ph = -ph;
    ph = ph * acc;
    if (bad == 1) ph = -ph;
    cs[2] = ph;
  }
}
// precision check pieces
__device__ inline bool isnanb(double x) { return (__double_as_longlong(x) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL && (__double_as_longlong(x) & 0x000FFFFFFFFFFFFFLL) != 0; }
__global__ void k_prec1(int n, const cd* A, const cd* B, double* pw, double* pm, double* pa, double* pb) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= n) return;
  double cs = 0, xm = 0, na = 0, nb = 0;
  for (int i = 0; i < n; i++) {
    cd a = A[i + (long)j * n], b = B[i + (long)j * n];
    bool fa = isnanb(a.r) || isnanb(a.i), fb = isnanb(b.r) || isnanb(b.i);
    if (fa) na += 1;
    if (fb) nb += 1;
    if (!(fa || fb)) { double d = cabsd(a - b); cs += d; xm = fmax(xm, d); }
  }
  pw[j] = cs; pm[j] = xm; pa[j] = na; pb[j] = nb;
}
__global__ void k_prec2(int n, const double* pw, const double* pm, const double* pa, const double* pb, double* r) {
  __shared__ double sh[1024];
  __shared__ double acc[3];
  if (threadIdx.x == 0) { acc[0] = 0; acc[1] = 0; acc[2] = 0; }
  __syncthreads();
  // ordered sum of pw through shared memory chunks
  double xs = 0;
  for (int c0 = 0; c0 < n; c0 += 1024) {
    int cn = min(1024, n - c0);
    for (int t = threadIdx.x; t < cn; t += blockDim.x) sh[t] = pw[c0 + t];
    __syncthreads();
    if (threadIdx.x == 0) { for (int t = 0; t < cn; t++) xs += sh[t]; }
    __syncthreads();
  }
  double xm = 0, na = 0, nb = 0;
  if (threadIdx.x == 0) { for (int t = 0; t < n; t++) { xm = fmax(xm, pm[t]); na += pa[t]; nb += pb[t]; }
    r[4] = xs; r[3] = xm; r[5] = na; r[6] = nb; }
}


// M(r,j) = (product over the factors of row r, in order) * M(r,j) -- as k_vslice_R of gpu_wrap_mod
__global__ void k_vslice_r(int nd, int ncol, int ld, int nf, int t, int cj, int desc, int ldptr, int nop, const int* vr_ptr,
                           const int* vr_ent, const int* V_meta, const int* V_type, const cd* nsig, const cd* V_exp, cd* M) {
  long tt = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (tt >= (long)nd * ncol) return;
  int r = tt % nd, j = tt / nd;
  int pa = vr_ptr[r + (long)(nf - 1) * ldptr], pb = vr_ptr[r + 1 + (long)(nf - 1) * ldptr] - 1;
  if (pb < pa) return;
  cd x = M[r + (long)j * ld];
  for (int q = 0; q <= pb - pa; q++) {
    int p = desc ? pb - q : pa + q;
    int n = vr_ent[2 * (p - 1)], i = vr_ent[2 * (p - 1) + 1];
    const int* mt = V_meta + 4L * ((n - 1) + (long)(nf - 1) * nop);
    if (mt[3] != 1) {
      int sp = (int)round(nsig[(n - 1) + (long)(t - 1) * nop].r);
      cd z = V_exp[mt[2] + i + (long)(sp + V_type[n - 1]) * mt[0] - 1];
      if (cj) z.i = -z.i;
      x = z * x;
    }
  }
  M[r + (long)j * ld] = x;
}

// ---- stabilisation helpers (gpu_stab_mod) and the column pre-pivoting of la_qrpp ------------------------------------
__global__ void k_ident(int n, int ld, cd* U) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; U[i + (long)j * ld] = (i == j) ? cd(1, 0) : cd(0, 0);
}
__global__ void k_setd(int n, cd* D) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) D[i] = cd(1, 0); }
__global__ void k_setl(int n, double* L) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) L[i] = 0.0; }
__global__ void k_diagabs(int n, int ld, const cd* U, cd* D) {
  int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
  cd u = U[i + (long)i * ld]; D[i] = cd(cabsd(u), 0.0);
}
// cs(1) = Phase = prod U(i,i) (* permutation sign, conjugated for side L), cs(2) = 1/Phase
__global__ void k_phase(int n, int ld, const cd* U, int sgn, int cj, cd* cs) {
  cd acc(1.0, 0.0);
  for (int i = threadIdx.x; i < n; i += 64) acc = acc * U[i + (long)i * ld];
  acc = wave_prod(acc);
  if (threadIdx.x == 0) {
    cd ph = acc; if (sgn < 0) ph = -ph; if (cj) ph = conjg(ph);
    cs[0] = ph; cs[1] = cdiv(cd(1.0, 0.0), ph);
  }
}
// cs(3) = conj(det A)/|det A| of the LU factors
__global__ void k_det_phase(int n, int ld, const cd* A, const int* ipiv, cd* cs) {
  cd acc(1.0, 0.0); int cnt = 0;
  for (int i = threadIdx.x; i < n; i += 64) { acc = acc * A[i + (long)i * ld]; if (ipiv[i] != i + 1) cnt++; }
  acc = wave_prod(acc);
  for (int o = 32; o > 0; o >>= 1) cnt += __shfl_xor(cnt, o, 64);
  if (threadIdx.x == 0) {
    cd det = acc; if (cnt & 1) det = -det; det = conjg(det);
    cs[2] = det / cabsd(det);
  }
}
__global__ void k_scale_row1(int n, int ld, cd* U, const cd* cs) { int j = blockIdx.x * blockDim.x + threadIdx.x; if (j < n) U[(long)j * ld] = cs[1] * U[(long)j * ld]; }
__global__ void k_scale_col1(int n, cd* U, const cd* cs) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) U[i] = cs[0] * U[i]; }
__global__ void k_triu(int n, int ld, const cd* U, int ldr, cd* R) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; R[i + (long)j * ldr] = (i <= j) ? U[i + (long)j * ld] : cd(0, 0);
}
// mode 0: W(i,j) = X(perm(i),j)  1: W(i,j) = X(i,perm(j))  2: W(perm(i),j) = X(i,j)  3: W(i,perm(j)) = X(i,j)
__global__ void k_perm(int mode, int n, const cd* X, int ldx, const int* perm, cd* W, int ldw) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n;
  if (mode == 0) W[i + (long)j * ldw] = X[(perm[i] - 1) + (long)j * ldx];
  else if (mode == 1) W[i + (long)j * ldw] = X[i + (long)(perm[j] - 1) * ldx];
  else if (mode == 2) W[(perm[i] - 1) + (long)j * ldw] = X[i + (long)j * ldx];
  else W[i + (long)(perm[j] - 1) * ldw] = X[i + (long)j * ldx];
}
__global__ void k_perm_vec(int n, const double* L, const int* perm, double* W) { int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) W[i] = L[perm[i] - 1]; }
__global__ void k_ct(int n, const cd* A, int lda, cd* B, int ldb) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; B[i + (long)j * ldb] = conjg(A[j + (long)i * lda]);
}
__global__ void k_colinv(int n, int lda, cd* A, const cd* D) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; double sv = 1.0 / D[j].r; A[i + (long)j * lda] = A[i + (long)j * lda] * sv;
}
// r(1) = number of NaN scales, r(2) = Re D(1), r(3) = Re D(n)
__global__ void k_testscale(int n, const cd* D, double* r) {
  __shared__ int s[256];
  int c = 0;
  for (int i = threadIdx.x; i < n; i += 256) if (isnanb(D[i].r) || isnanb(D[i].i)) c++;
  s[threadIdx.x] = c; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) s[threadIdx.x] += s[threadIdx.x + o]; __syncthreads(); }
  if (threadIdx.x == 0) { r[0] = (double)s[0]; r[1] = D[0].r; r[2] = D[n - 1].r; }
}
// nrm(j) = 2-norm of column j (scaled); NaN if the column holds a NaN or Inf; one block per column
__global__ void k_colnorms(int m, int n, const cd* A, int lda, double* nrm) {
  __shared__ double sm[256]; __shared__ int sf[256];
  int j = blockIdx.x; const cd* a = A + (long)j * lda;
  double mx = 0; int nf = 0;
  for (int i = threadIdx.x; i < m; i += 256) {
    double re = a[i].r, im = a[i].i;
    if (((__double_as_longlong(re) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL) || ((__double_as_longlong(im) & 0x7FF0000000000000LL) == 0x7FF0000000000000LL)) nf++;
    mx = fmax(mx, fmax(fabs(re), fabs(im)));
  }
  sm[threadIdx.x] = mx; sf[threadIdx.x] = nf; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) { sm[threadIdx.x] = fmax(sm[threadIdx.x], sm[threadIdx.x + o]); sf[threadIdx.x] += sf[threadIdx.x + o]; } __syncthreads(); }
  mx = sm[0]; nf = sf[0]; __syncthreads();
  double s = 0;
  if (nf == 0 && mx > 0.0) for (int i = threadIdx.x; i < m; i += 256) { double v = a[i].r / mx; s += v * v; v = a[i].i / mx; s += v * v; }
  sm[threadIdx.x] = s; __syncthreads();
  for (int o = 128; o > 0; o >>= 1) { if (threadIdx.x < o) sm[threadIdx.x] += sm[threadIdx.x + o]; __syncthreads(); }
  if (threadIdx.x == 0) nrm[j] = (nf > 0) ? __longlong_as_double(0x7FF8000000000000LL) : (mx > 0.0 ? mx * sqrt(sm[0]) : 0.0);
}
__global__ void k_gather_cols(int m, int n, const cd* A, int lda, const int* ipvt, cd* scr) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)m * n) return;
  int i = t % m, j = t / m; scr[t] = A[i + (long)(ipvt[j] - 1) * lda];
}
__global__ void k_scatter_cols(int m, int n, const cd* scr, cd* A, int lda) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)m * n) return;
  int i = t % m, j = t / m; A[i + (long)j * lda] = scr[t];
}

extern "C" {

int alfhk_ident(int n, int ld, void* U) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_ident, nb(t,256), 256, n, ld, (cd*)U); return fin(); }
int alfhk_setd(int n, void* D) { if (n<=0) return 0; LAUNCH(k_setd, nb(n,256), 256, n, (cd*)D); return fin(); }
int alfhk_setl(int n, void* L) { if (n<=0) return 0; LAUNCH(k_setl, nb(n,256), 256, n, (double*)L); return fin(); }
int alfhk_diagabs(int n, int ld, void* U, void* D) { if (n<=0) return 0; LAUNCH(k_diagabs, nb(n,256), 256, n, ld, (const cd*)U, (cd*)D); return fin(); }
int alfhk_phase(int n, int ld, void* U, int sgn, int cj, void* cs) { LAUNCH(k_phase, 1, 64, n, ld, (const cd*)U, sgn, cj, (cd*)cs); return fin(); }
int alfhk_det_phase(int n, int ld, void* A, void* ipiv, void* cs) { LAUNCH(k_det_phase, 1, 64, n, ld, (const cd*)A, (const int*)ipiv, (cd*)cs); return fin(); }
int alfhk_scale_row1(int n, int ld, void* U, void* cs) { if (n<=0) return 0; LAUNCH(k_scale_row1, nb(n,256), 256, n, ld, (cd*)U, (const cd*)cs); return fin(); }
int alfhk_scale_col1(int n, void* U, void* cs) { if (n<=0) return 0; LAUNCH(k_scale_col1, nb(n,256), 256, n, (cd*)U, (const cd*)cs); return fin(); }
int alfhk_triu(int n, int ld, void* U, int ldr, void* R) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_triu, nb(t,256), 256, n, ld, (const cd*)U, ldr, (cd*)R); return fin(); }
int alfhk_perm(int mode, int n, void* X, int ldx, void* perm, void* W, int ldw) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_perm, nb(t,256), 256, mode, n, (const cd*)X, ldx, (const int*)perm, (cd*)W, ldw); return fin(); }
int alfhk_perm_vec(int n, void* L, void* perm, void* W) { if (n<=0) return 0; LAUNCH(k_perm_vec, nb(n,256), 256, n, (const double*)L, (const int*)perm, (double*)W); return fin(); }
int alfhk_ct(int n, void* A, int lda, void* B, int ldb) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_ct, nb(t,256), 256, n, (const cd*)A, lda, (cd*)B, ldb); return fin(); }
int alfhk_colinv(int n, int lda, void* A, void* D) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_colinv, nb(t,256), 256, n, lda, (cd*)A, (const cd*)D); return fin(); }
int alfhk_testscale(int n, void* D, void* r) { if (n<=0) return 0; LAUNCH(k_testscale, 1, 256, n, (const cd*)D, (double*)r); return fin(); }
int alfhk_colnorms(int m, int n, void* A, int lda, void* nrm) { if (n<=0) return 0; LAUNCH(k_colnorms, n, 256, m, n, (const cd*)A, lda, (double*)nrm); return fin(); }
static cd* gscr = nullptr; static long gscr_n = 0;
int alfhk_permute_cols(int m, int n, void* A, int lda, void* ipvt) {
  long t=(long)m*n; if (t<=0) return 0;
  if (gscr_n < t) { if (gscr) hipFree(gscr); gscr = nullptr; gscr_n = 0; if (hipMalloc((void**)&gscr, t * sizeof(cd)) != hipSuccess) return 3; gscr_n = t; }
  LAUNCH(k_gather_cols, nb(t,256), 256, m, n, (const cd*)A, lda, (const int*)ipvt, gscr);
  LAUNCH(k_scatter_cols, nb(t,256), 256, m, n, (const cd*)gscr, (cd*)A, lda);
  return fin();
}

int alfhk_copyz(long n, void* a, void* b) {
  if (n <= 0) return 0;
  LAUNCH(k_copyz, nb(n, 256), 256, n, (const cd*)a, (cd*)b);
  return fin();
}
int alfhk_blk_gather(int nfl, int nfa, int nfb, int b, int nd, int ld, void* G, void* bsd, void* gbd) {
  long tot = (long)nfl * b * b; if (tot <= 0) return 0;
  LAUNCH(k_blk_gather, nb(tot, 128), 128, nfl, nfa, nfb, b, nd, ld, (const cd*)G, (const int*)bsd, (cd*)gbd);
  return fin();
}
int alfhk_blk_scale_rows(int nf, int fl, int ld, int b, int nd, int cbld, void* G, void* bsd, void* cbuf) {
  long tot = (long)b * ld; if (tot <= 0) return 0;
  LAUNCH(k_blk_scale_rows, nb(tot, 256), 256, nf, fl, ld, b, nd, cbld, (cd*)G, (const int*)bsd, (const cd*)cbuf);
  return fin();
}
int alfhk_blk_scale_cols(int nf, int fl, int ld, int b, int nd, int cbld, void* G, void* bsd, void* cbuf) {
  long tot = (long)b * ld; if (tot <= 0) return 0;
  LAUNCH(k_blk_scale_cols, nb(tot, 256), 256, nf, fl, ld, b, nd, cbld, (cd*)G, (const int*)bsd, (const cd*)cbuf);
  return fin();
}
int alfhk_blk_end(int nf, int fl, int ld, int b, int nd, int cbld, int kk, int nslot_ld, void* G, void* bsd, void* cbuf, void* S, void* C, void* Qt, void* slot) {
  long tot = 2L * b * ld; if (tot <= 0) return 0;
  LAUNCH(k_blk_end, nb(tot, 256), 256, nf, fl, ld, b, nd, cbld, kk, nslot_ld, (cd*)G, (const int*)bsd, (const cd*)cbuf, (const int*)S, (cd*)C, (cd*)Qt, (int*)slot);
  return fin();
}
int alfhk_flush_gather(int ld, int nrow, int k, int nf, int nslot_ld, void* G, void* S, void* C, void* Qt, void* slot) {
  long tot = (long)k * nrow; if (tot <= 0) return 0;
  LAUNCH(k_flush_gather, nb(tot, 256), 256, ld, nrow, k, nf, nslot_ld, (const cd*)G, (const int*)S, (cd*)C, (cd*)Qt, (int*)slot);
  return fin();
}
int alfhk_scatter_rows(int nd, int ncol, void* P, void* M, int ld, void* W, int ldw) {
  long tot = (long)nd * ncol; if (tot <= 0) return 0;
  LAUNCH(k_scatter_rows, nb(tot, 256), 256, nd, ncol, (const int*)P, (cd*)M, ld, (const cd*)W, ldw);
  return fin();
}
int alfhk_corr(int kind, int k, int N, int Norb, int NFL, double ns, double sq3, double wr, double wi, int ld,
               void* G, void* jc, void* inv, void* acc, void* acc0) {
  long MT = (long)N * Norb * Norb;
  long blocks = MT + ((kind == 3 || kind == 4) ? Norb : 0);
  LAUNCH(k_corr, blocks, RB, kind, k, N, Norb, NFL, ns, sq3, wr, wi, ld, (const cd*)G, (const int*)jc, (const int*)inv, (cd*)acc, (cd*)acc0);
  return fin();
}
static double* dres = nullptr;
int alfhk_scalars(int nk, void* hi, void* hj, void* hf, void* hc, int N, int Norb, int NFL, int dec, int ld, void* G, void* inv,
                  void* uo, double* hres) {
  if (!dres && hipMalloc((void**)&dres, 6 * sizeof(double)) != hipSuccess) return 3;
  LAUNCH(k_scalars, 1, RB, nk, (const int*)hi, (const int*)hj, (const int*)hf, (const cd*)hc, N, Norb, NFL, dec, ld, (const cd*)G,
         (const int*)inv, (const double*)uo, dres);
  int r = fin();
  if (r) return r;
  return hipMemcpy(hres, dres, 6 * sizeof(double), hipMemcpyDeviceToHost) == hipSuccess ? 0 : 4;
}
int alfhk_local(int N, int Norb, int lNt, int it, int dec, double wt, double wr, double wi, int ld, void* inv, void* G, void* lacc1,
                void* lacc2) {
  long tot = (long)N * Norb; if (tot <= 0) return 0;
  LAUNCH(k_local, nb(tot, 256), 256, N, Norb, lNt, it, dec, wt, wr, wi, ld, (const int*)inv, (const cd*)G, (cd*)lacc1, (cd*)lacc2);
  return fin();
}


int alfhk_divrows(int n, int ld, void* U, void* D) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_divrows, nb(t,256), 256, n, ld, (cd*)U, (const cd*)D); return fin(); }
int alfhk_rowdiv(int n, int lda, void* A, void* D) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_rowdiv, nb(t,256), 256, n, lda, (cd*)A, (const cd*)D); return fin(); }
int alfhk_scale_cols(int n, int ld, void* U, void* D) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_scale_cols, nb(t,256), 256, n, ld, (cd*)U, (const cd*)D); return fin(); }
int alfhk_tpup(int n, void* T, void* R, void* Dr, void* Dl) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_tpup, nb(t,256), 256, n, (cd*)T, (const cd*)R, (const cd*)Dr, (const cd*)Dl); return fin(); }
int alfhk_ns(int n, void* A) { long t=(long)n*n; if (t<=0) return 0;
  LAUNCH(k_ns, nb(t,256), 256, n, (cd*)A); return fin(); }
int alfhk_cgr_phase(int n, int ld, void* T, void* tau, void* cs, int nvar, int sgn, int bad) {
  LAUNCH(k_cgr_phase, 1, 64, n, ld, (const cd*)T, (const cd*)tau, (cd*)cs, nvar, sgn, bad); return fin(); }
static double* pscr = nullptr; static long pscr_n = 0;
int alfhk_prec(int n, void* A, void* B, void* r) {
  if (pscr_n < n) { if (pscr) hipFree(pscr); if (hipMalloc((void**)&pscr, 4L * n * sizeof(double)) != hipSuccess) return 3; pscr_n = n; }
  double *pw = pscr, *pm = pscr + n, *pa = pscr + 2L * n, *pb = pscr + 3L * n;
  LAUNCH(k_prec1, nb(n, 64), 64, n, (const cd*)A, (const cd*)B, pw, pm, pa, pb);
  LAUNCH(k_prec2, 1, 256, n, (const double*)pw, (const double*)pm, (const double*)pa, (const double*)pb, (double*)r);
  return fin();
}

int alfhk_vslice_r(int nd, int ncol, int ld, int nf, int t, int cj, int desc, int ldptr, int nop, void* vr_ptr, void* vr_ent,
                   void* V_meta, void* V_type, void* nsig, void* V_exp, void* M) {
  long tot = (long)nd * ncol; if (tot <= 0) return 0;
  LAUNCH(k_vslice_r, nb(tot, 256), 256, nd, ncol, ld, nf, t, cj, desc, ldptr, nop, (const int*)vr_ptr, (const int*)vr_ent,
         (const int*)V_meta, (const int*)V_type, (const cd*)nsig, (const cd*)V_exp, (cd*)M);
  return fin();
}
}  // extern "C"
