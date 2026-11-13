// GPU=apple (ALF_LA_METAL): Apple-GPU assisted ZGEMM for the CPU build.
//
// Apple GPUs have no float64. Large ZGEMMs (all of M, N, K >= ALF_METAL_MIN, default 512) are split by
// columns of C: the GPU computes the first ~ALF_METAL_FRAC (default 0.35) of the columns in float-float
// ("double-single", ~48-bit) arithmetic while Accelerate computes the rest on the CPU at the same time.
// Smaller or unsupported calls go to Accelerate unchanged. Measured (M1 Pro, plan/APPLE_GPU_FEASIBILITY.md):
// normwise error ~1e-12 relative to max|C| (vs ~1e-16 in double), 1.2-1.45x faster than Accelerate alone at N >= 576.
//
// This file defines the Fortran symbol zgemm_, so it replaces the BLAS ZGEMM for all ALF Fortran code
// (calls made inside LAPACK are not affected). Inputs outside [2^-40, 2^40] in max-norm fall back to the CPU.
// Env: ALF_METAL_MIN (int), ALF_METAL_FRAC (0..0.9; 0 disables the GPU), ALF_METAL_VERBOSE=1 (stats at exit),
//      ALF_METAL_FASTMATH=1 (planted-defect control: compiles the kernels with fast math; never for production).
#import <Metal/Metal.h>
// Accelerate's legacy CBLAS entry point, declared by hand: including the Accelerate headers would also declare
// zgemm_ (with ACCELERATE_NEW_LAPACK even under another symbol name), which this file defines.
enum { CblasColMajor = 102, CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113 };
extern void cblas_zgemm(int order, int ta, int tb, int m, int n, int k, const void* alpha, const void* A, int lda,
                        const void* B, int ldb, const void* beta, void* C, int ldc);
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <dispatch/dispatch.h>

typedef double _Complex zc;

static const char* ksrc =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"#define RM 2\n#define RN 2\n#define TX 8\n#define TY 8\n#define BK 16\n"
"static inline float2 d2ff(uint2 u){ uint hi32=u.y, lo32=u.x; int e=int((hi32>>20)&0x7FFu); if(e==0) return float2(0.0f);\n"
"  float s=(hi32>>31)?-1.0f:1.0f; uint m23=((hi32&0xFFFFFu)<<3)|(lo32>>29); float h=as_type<float>(0x3F800000u|m23);\n"
"  float l=float(lo32&0x1FFFFFFFu)*0x1p-52f; int ex=e-1023; if(ex>100) return float2(INFINITY); if(ex<-120) return float2(0.0f);\n"
"  return float2(s*ldexp(h,ex), s*ldexp(l,ex)); }\n"
"struct ConvP { uint rows, cols, ld, trans, conj, off; };\n"
// dst (rows x cols, packed col-major) = op(X); X raw doubles (re,im pairs) at complex-element offset off, leading dim ld
"kernel void conv_op(device const uint2* src [[buffer(0)]], device float4* dst [[buffer(1)]], constant ConvP& p [[buffer(2)]],\n"
"                    device atomic_uint* mx [[buffer(3)]], uint2 g [[thread_position_in_grid]]){\n"
"  uint i=g.x, j=g.y; if(i>=p.rows||j>=p.cols) return;\n"
"  ulong idx = p.trans ? (ulong)i*p.ld + j : (ulong)j*p.ld + i; idx += p.off;\n"
"  float2 re=d2ff(src[2*idx]), im=d2ff(src[2*idx+1]); if(p.conj){ im=-im; }\n"
"  dst[(ulong)j*p.rows+i]=float4(re.x,re.y,im.x,im.y);\n"
"  float a=max(abs(re.x),abs(im.x)); if(a>0) atomic_fetch_max_explicit(mx, as_type<uint>(a), memory_order_relaxed); }\n"
"static inline void ffmac(thread float2& s, float2 a, float2 b){ float p=a.x*b.x; float e=fma(a.x,b.x,-p); e=fma(a.x,b.y,e); e=fma(a.y,b.x,e);\n"
"  float t=s.x+p; float bb=t-s.x; float err=(s.x-(t-bb))+(p-bb); s.x=t; s.y+=err+e; }\n"
"#define BM (TX*RM)\n#define BN (TY*RN)\n"
"kernel void zgemm_ff(device const float4* A [[buffer(0)]], device const float4* B [[buffer(1)]], device float4* C [[buffer(2)]],\n"
"  constant uint4& dims [[buffer(3)]], uint2 tg [[threadgroup_position_in_grid]], uint2 tid [[thread_position_in_threadgroup]]){\n"
"  const uint M=dims.x, N=dims.y, K=dims.z; threadgroup float4 As[BK][BM]; threadgroup float4 Bs[BK][BN];\n"
"  const uint tx=tid.x, ty=tid.y, lin=ty*TX+tx; const uint m0=tg.x*BM, n0=tg.y*BN; float2 cr[RM][RN], ci[RM][RN];\n"
"  for(int i=0;i<RM;i++) for(int j=0;j<RN;j++){ cr[i][j]=float2(0); ci[i][j]=float2(0); }\n"
"  for(uint k0=0;k0<K;k0+=BK){\n"
"    for(uint idx=lin; idx<BK*BM; idx+=TX*TY){ uint r=idx%BM, kk=idx/BM; uint gm=m0+r, gk=k0+kk; As[kk][r]=(gm<M&&gk<K)?A[(ulong)gk*M+gm]:float4(0); }\n"
"    for(uint idx=lin; idx<BK*BN; idx+=TX*TY){ uint kk=idx%BK, c=idx/BK; uint gk=k0+kk, gn=n0+c; Bs[kk][c]=(gk<K&&gn<N)?B[(ulong)gn*K+gk]:float4(0); }\n"
"    threadgroup_barrier(mem_flags::mem_threadgroup);\n"
"    for(int kk=0;kk<BK;kk++){ float4 a[RM], b[RN]; for(int i=0;i<RM;i++) a[i]=As[kk][tx+i*TX]; for(int j=0;j<RN;j++) b[j]=Bs[kk][ty+j*TY];\n"
"      for(int i=0;i<RM;i++) for(int j=0;j<RN;j++){ ffmac(cr[i][j],a[i].xy,b[j].xy); ffmac(cr[i][j],a[i].zw,-b[j].zw); ffmac(ci[i][j],a[i].xy,b[j].zw); ffmac(ci[i][j],a[i].zw,b[j].xy); } }\n"
"    threadgroup_barrier(mem_flags::mem_threadgroup); }\n"
"  for(int i=0;i<RM;i++) for(int j=0;j<RN;j++){ uint gm=m0+tx+i*TX, gn=n0+ty+j*TY; if(gm<M&&gn<N) C[(ulong)gn*M+gm]=float4(cr[i][j].x,cr[i][j].y,ci[i][j].x,ci[i][j].y); } }\n";

struct ConvP { uint32_t rows, cols, ld, trans, conj, off; };

static int g_init = 0, g_ok = 0, g_min = 512, g_verbose = 0;
static double g_frac = 0.35;
static id<MTLDevice> g_dev; static id<MTLCommandQueue> g_q; static id<MTLComputePipelineState> g_pconv, g_pgemm;
static id<MTLBuffer> g_wA, g_wB, g_wC, g_mx; static size_t g_cap = 0;
static long n_gpu = 0, n_cpu = 0, n_fallback = 0, n_big = 0; static double t_gpu = 0, t_cpu = 0, fl_gpu = 0, t_big = 0, fl_big = 0;

static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void report(void){
  if (!g_verbose) return;
  fprintf(stderr, "ALF_METAL zgemm: %ld split calls (%.3f s, %.1f GF/s, %ld CPU fallbacks), %ld CPU-only calls (%.3f s;"
          " of these %ld above the size threshold: %.3f s, %.1f GF/s)\n",
          n_gpu, t_gpu, t_gpu > 0 ? fl_gpu / t_gpu * 1e-9 : 0.0, n_fallback, n_cpu, t_cpu, n_big, t_big, t_big > 0 ? fl_big / t_big * 1e-9 : 0.0);
}
static void init(void){
  g_init = 1;
  const char* e;
  if ((e = getenv("ALF_METAL_MIN"))) g_min = atoi(e);
  if ((e = getenv("ALF_METAL_FRAC"))) g_frac = atof(e);
  if ((e = getenv("ALF_METAL_VERBOSE"))) g_verbose = atoi(e);
  atexit(report);
  if (g_frac <= 0) return;
  @autoreleasepool {
    g_dev = MTLCreateSystemDefaultDevice(); if (!g_dev) return;
    g_q = [g_dev newCommandQueue];
    MTLCompileOptions* o = [MTLCompileOptions new];
    // float-float arithmetic needs IEEE semantics: no fast math, no reassociation
    // ALF_METAL_FASTMATH=1 is a planted-defect control only (fast math destroys float-float: errors ~1e-6)
    int fast = getenv("ALF_METAL_FASTMATH") && atoi(getenv("ALF_METAL_FASTMATH"));
    if (@available(macOS 15.0, *)) { o.mathMode = fast ? MTLMathModeFast : MTLMathModeSafe; o.mathFloatingPointFunctions = fast ? MTLMathFloatingPointFunctionsFast : MTLMathFloatingPointFunctionsPrecise; }
    NSError* err = nil;
    id<MTLLibrary> l = [g_dev newLibraryWithSource:[NSString stringWithUTF8String:ksrc] options:o error:&err];
    if (!l) { fprintf(stderr, "ALF_METAL: shader compile failed: %s\n", [[err description] UTF8String]); return; }
    g_pconv = [g_dev newComputePipelineStateWithFunction:[l newFunctionWithName:@"conv_op"] error:&err];
    g_pgemm = [g_dev newComputePipelineStateWithFunction:[l newFunctionWithName:@"zgemm_ff"] error:&err];
    if (!g_pconv || !g_pgemm) return;
    g_mx = [g_dev newBufferWithLength:16 options:MTLResourceStorageModeShared];
    g_ok = 1;
    if (g_verbose) fprintf(stderr, "ALF_METAL: %s, min=%d frac=%.2f\n", [[g_dev name] UTF8String], g_min, g_frac);
  }
}
static void ensure(size_t a, size_t b, size_t c){
  size_t m = a > b ? a : b; m = m > c ? m : c;
  if (m <= g_cap) return;
  g_cap = m;
  g_wA = [g_dev newBufferWithLength:16 * m options:MTLResourceStorageModeShared];
  g_wB = [g_dev newBufferWithLength:16 * m options:MTLResourceStorageModeShared];
  g_wC = [g_dev newBufferWithLength:16 * m options:MTLResourceStorageModeShared];
}
// Wrap host memory [p, p+bytes) without copy: page-align the start down, return the offset in complex elements.
// Returns nil (caller copies) if p is not 16-byte aligned.
static id<MTLBuffer> wrap(const void* p, size_t bytes, uint32_t* off_z){
  uintptr_t pg = 16384, a = (uintptr_t)p, base = a & ~(pg - 1), end = (a + bytes + pg - 1) & ~(pg - 1);
  if (a & 15) return nil;
  *off_z = (uint32_t)((a - base) / 16);
  id<MTLBuffer> b = [g_dev newBufferWithBytesNoCopy:(void*)base length:end - base options:MTLResourceStorageModeShared deallocator:nil];
  return b;
}
static void cpu_zgemm(char ta, char tb, int m, int n, int k, const zc* al, const zc* A, int lda, const zc* B, int ldb, const zc* be, zc* C, int ldc){
  int TA = (ta == 'N' || ta == 'n') ? CblasNoTrans : (ta == 'T' || ta == 't') ? CblasTrans : CblasConjTrans;
  int TB = (tb == 'N' || tb == 'n') ? CblasNoTrans : (tb == 'T' || tb == 't') ? CblasTrans : CblasConjTrans;
  cblas_zgemm(CblasColMajor, TA, TB, m, n, k, al, A, lda, B, ldb, be, C, ldc);
}

void zgemm_(const char* ta, const char* tb, const int* M, const int* N, const int* K, const zc* alpha,
            const zc* A, const int* lda, const zc* B, const int* ldb, const zc* beta, zc* C, const int* ldc,
            size_t lta, size_t ltb){
  (void)lta; (void)ltb;
  if (!g_init) init();  // not thread-safe: ALF calls ZGEMM from one thread
  int m = *M, n = *N, k = *K;
  double t0 = g_verbose ? now() : 0;
  if (!g_ok || m < g_min || n < g_min || k < g_min) {
    cpu_zgemm(*ta, *tb, m, n, k, alpha, A, *lda, B, *ldb, beta, C, *ldc);
    if (g_verbose) { double dt = now() - t0; n_cpu++; t_cpu += dt;
      if (m >= g_min && n >= g_min && k >= g_min) { n_big++; t_big += dt; fl_big += 8.0 * m * n * k; } }
    return;
  }
  int ng = ((int)(n * g_frac)) & ~31; if (ng < 32) ng = 32;
  int tA = !(*ta == 'N' || *ta == 'n'), cA = (*ta == 'C' || *ta == 'c');
  int tB = !(*tb == 'N' || *tb == 'n'), cB = (*tb == 'C' || *tb == 'c');
  int LA = *lda, LB = *ldb, LC = *ldc;
  zc al = *alpha, be = *beta;
  @autoreleasepool {
    ensure((size_t)m * k, (size_t)k * ng, (size_t)m * ng);
    // op(A): m x k; storage A is (tA ? k x m : m x k); op(B)(:,0:ng): k x ng
    size_t bytesA = 16 * ((size_t)LA * ((tA ? m : k) - 1) + (tA ? k : m));
    size_t bytesB = 16 * ((size_t)LB * ((tB ? k : ng) - 1) + (tB ? ng : k));
    struct ConvP pa = { (uint32_t)m, (uint32_t)k, (uint32_t)LA, (uint32_t)tA, (uint32_t)cA, 0 };
    struct ConvP pb = { (uint32_t)k, (uint32_t)ng, (uint32_t)LB, (uint32_t)tB, (uint32_t)cB, 0 };
    id<MTLBuffer> bA = wrap(A, bytesA, &pa.off), bB = wrap(B, bytesB, &pb.off);
    if (!bA) bA = [g_dev newBufferWithBytes:A length:bytesA options:MTLResourceStorageModeShared], pa.off = 0;
    if (!bB) bB = [g_dev newBufferWithBytes:B length:bytesB options:MTLResourceStorageModeShared], pb.off = 0;
    uint32_t* mx = (uint32_t*)g_mx.contents; mx[0] = 0; mx[1] = 0;
    id<MTLCommandBuffer> cb = [g_q commandBuffer];
    id<MTLComputeCommandEncoder> en = [cb computeCommandEncoder];
    [en setComputePipelineState:g_pconv];
    [en setBuffer:bA offset:0 atIndex:0]; [en setBuffer:g_wA offset:0 atIndex:1]; [en setBytes:&pa length:sizeof pa atIndex:2]; [en setBuffer:g_mx offset:0 atIndex:3];
    [en dispatchThreads:MTLSizeMake(m, k, 1) threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
    [en setBuffer:bB offset:0 atIndex:0]; [en setBuffer:g_wB offset:0 atIndex:1]; [en setBytes:&pb length:sizeof pb atIndex:2]; [en setBuffer:g_mx offset:4 atIndex:3];
    [en dispatchThreads:MTLSizeMake(k, ng, 1) threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
    uint32_t dims[4] = { (uint32_t)m, (uint32_t)ng, (uint32_t)k, 0 };
    [en setComputePipelineState:g_pgemm];
    [en setBuffer:g_wA offset:0 atIndex:0]; [en setBuffer:g_wB offset:0 atIndex:1]; [en setBuffer:g_wC offset:0 atIndex:2]; [en setBytes:dims length:16 atIndex:3];
    [en dispatchThreadgroups:MTLSizeMake((m + 15) / 16, (ng + 15) / 16, 1) threadsPerThreadgroup:MTLSizeMake(8, 8, 1)];
    [en endEncoding];
    [cb commit];
    // CPU part, concurrently: columns ng..n-1
    if (ng < n) {
      const zc* Bs = tB ? B + ng : B + (size_t)LB * ng;
      cpu_zgemm(*ta, *tb, m, n - ng, k, alpha, A, LA, Bs, LB, beta, C + (size_t)LC * ng, LC);
    }
    [cb waitUntilCompleted];
    float fa, fb; memcpy(&fa, &mx[0], 4); memcpy(&fb, &mx[1], 4);
    int okrange = (fa >= 0x1p-40f && fa <= 0x1p40f) && (fb >= 0x1p-40f && fb <= 0x1p40f) && cb.status == MTLCommandBufferStatusCompleted;
    if (!okrange) {  // values the float-float path cannot represent well (or a zero matrix): redo on the CPU
      cpu_zgemm(*ta, *tb, m, ng, k, alpha, A, LA, B, LB, beta, C, LC);
      n_fallback++;
    } else {
      const float* P = (const float*)g_wC.contents;
      int nchunk = 8;
      dispatch_apply(nchunk, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t c){
        int j0 = (int)(ng * c / nchunk), j1 = (int)(ng * (c + 1) / nchunk);
        for (int j = j0; j < j1; j++) {
          zc* Cj = C + (size_t)LC * j; const float* Pj = P + 4 * (size_t)m * j;
          if (be == 0) for (int i = 0; i < m; i++) Cj[i] = al * (((double)Pj[4*i] + (double)Pj[4*i+1]) + I * ((double)Pj[4*i+2] + (double)Pj[4*i+3]));
          else         for (int i = 0; i < m; i++) Cj[i] = al * (((double)Pj[4*i] + (double)Pj[4*i+1]) + I * ((double)Pj[4*i+2] + (double)Pj[4*i+3])) + be * Cj[i];
        }
      });
    }
  }
  if (g_verbose) { n_gpu++; t_gpu += now() - t0; fl_gpu += 8.0 * m * n * k; }
}
