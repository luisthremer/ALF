import ctypes as C, sys, h5py, numpy as np
ENV="/home/atuin/b133ae/b133ae29/SOFTWARE/micromamba/envs/general_python/lib/"
ob=C.CDLL(ENV+"libblas.so.3"); lm=C.CDLL("libm.so.6")
ob.cblas_dgemm.restype=None; ob.cblas_zgemm.restype=None
class Z(C.Structure): _fields_=[("re",C.c_double),("im",C.c_double)]
lm.cexp.restype=Z; lm.cexp.argtypes=[Z]
P=lambda a: a.ctypes.data_as(C.c_void_p)
I=C.c_int
def dgemm(order,ta,tb,M,N,K,A,lda,B,ldb,ldc):
    Cm=np.zeros(M*N); ob.cblas_dgemm(I(order),I(ta),I(tb),I(M),I(N),I(K),C.c_double(1.0),P(A),I(lda),P(B),I(ldb),C.c_double(0.0),P(Cm),I(ldc)); return Cm
def zgemm(order,ta,tb,M,N,K,A,lda,B,ldb,ldc):
    Cm=np.zeros(M*N,complex); one=np.array([1+0j]); zero=np.array([0j])
    ob.cblas_zgemm(I(order),I(ta),I(tb),I(M),I(N),I(K),P(one),P(A),I(lda),P(B),I(ldb),P(zero),P(Cm),I(ldc)); return Cm
RM,CM,NT,T=101,102,111,112
ok=True
for f in sys.argv[1:]:
    g=h5py.File(f,"r")["SpinZ_localtau"]
    spin=g["obser"][:3,0,:,:,0]; L=int(np.sqrt(spin.shape[2])); Ns=L*L; Nt=spin.shape[1]
    def ll(i):
        o=L//2-1; return (i//L)-o,(i%L)-o
    n=np.arange(-L//2+1,L//2+1); q=2*np.pi/L*n
    rl=np.array([(a,b) for a in q for b in q]); lt=np.array([ll(i) for i in range(Ns)])
    qr=rl@lt.T; ph=np.exp(-1j*qr); ref=spin@ph.T
    ltf=np.ascontiguousarray(lt,dtype=float)        # (Ns,2) row-major
    qa=dgemm(RM,NT,T,Ns,Ns,2,rl,2,ltf,2,Ns).reshape(Ns,Ns)
    ltT=np.ascontiguousarray(lt.T,dtype=float)      # (2,Ns) row-major
    qb=dgemm(RM,NT,NT,Ns,Ns,2,rl,2,ltT,Ns,Ns).reshape(Ns,Ns)
    phc=np.array([[ (lambda z: complex(z.re,z.im))(lm.cexp(Z(0.0,-x))) for x in row] for row in qr])
    sq=np.empty_like(ref)
    for b in range(spin.shape[0]):
        A=np.ascontiguousarray(spin[b]).astype(complex)
        sq[b]=zgemm(RM,NT,T,Nt,Ns,Ns,A,Ns,ph,Ns,Ns).reshape(Nt,Ns)
    r=(L,np.array_equal(qa,qr),np.array_equal(qb,qr),np.array_equal(phc,ph),np.array_equal(sq,ref))
    print("L=%d qr(RM,N,T)=%s qr(RM,N,N)=%s cexp=%s zgemm(RM,N,T)=%s"%r); ok&=all(r[1:2]) and r[3] and r[4]
print("ALL OK" if ok else "MISMATCH")
