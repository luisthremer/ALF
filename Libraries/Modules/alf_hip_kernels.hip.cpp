// ALF GPU=amd: small device kernels written in HIP instead of OpenMP target regions.
// A libgomp-GCN target region with a parallel loop costs ~150 us to launch; a HIP kernel launch + hipDeviceSynchronize ~15 us.
// Every function takes device pointers (from use_device_addr in alf_hipk_mod / linalg_hip_mod) and ends in fin() = alf_hip_fin(): the work is synchronised
// before it returns, or, with the deferred synchronisation (below), only marked as pending and synchronised by alf_hip_presync() before
// the next OpenMP kernel launch or device copy, because the OpenMP target regions run on their own queues. Index conventions follow
// the Fortran arrays (column major).
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


// ---- linalg_hip_mod helpers: copies, identity, column pre-pivoting of la_zqrpp --------------------------------
__global__ void k_copyz(long n, const cd* a, cd* b) {
  long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) b[i] = a[i];
}
__global__ void k_ident(int n, int ld, cd* U) {
  long t = (long)blockIdx.x * blockDim.x + threadIdx.x; if (t >= (long)n * n) return;
  int i = t % n, j = t / n; U[i + (long)j * ld] = (i == j) ? cd(1, 0) : cd(0, 0);
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
int alfhk_copyz(long n, void* a, void* b) {
  if (n <= 0) return 0;
  LAUNCH(k_copyz, nb(n, 256), 256, n, (const cd*)a, (cd*)b);
  return fin();
}
int alfhk_ident(int n, int ld, void* U) { long t=(long)n*n; if (t<=0) return 0; LAUNCH(k_ident, nb(t,256), 256, n, ld, (cd*)U); return fin(); }
int alfhk_colnorms(int m, int n, void* A, int lda, void* nrm) { if (n<=0) return 0; LAUNCH(k_colnorms, n, 256, m, n, (const cd*)A, lda, (double*)nrm); return fin(); }
static cd* gscr = nullptr; static long gscr_n = 0;
int alfhk_permute_cols(int m, int n, void* A, int lda, void* ipvt) {
  long t=(long)m*n; if (t<=0) return 0;
  if (gscr_n < t) { if (gscr) (void)hipFree(gscr); gscr = nullptr; gscr_n = 0; if (hipMalloc((void**)&gscr, t * sizeof(cd)) != hipSuccess) return 3; gscr_n = t; }
  LAUNCH(k_gather_cols, nb(t,256), 256, m, n, (const cd*)A, lda, (const int*)ipvt, gscr);
  LAUNCH(k_scatter_cols, nb(t,256), 256, m, n, (const cd*)gscr, (cd*)A, lda);
  return fin();
}
}  // extern "C"
