import numpy as np, re, subprocess, sys
from scipy.linalg import expm
N=16
O=np.zeros((N,N),complex)
for l in open(sys.argv[1]):
    if l.startswith("DBGO"):
        f=l.split(); O[int(f[1])-1,int(f[2])-1]=float(f[3])+1j*float(f[4])
dt=0.1; L=40; pin=5   # site 6 (1-based)
E=np.zeros((N,N)); E[pin,pin]=-1.0     # flavour 1: O2 = -h0 on pin, h0=1 ; g_t=-0.5 multiplies exp(g*O2)
def run(shift, symm):
    g=np.zeros(L+2); g[20+shift]=-0.5; g[21+shift]=-0.5
    def slot(t):
        gt = 0.5*(g[t]+(g[t-1] if t>1 else g[L])) if symm else g[t]
        return expm(-dt*O)@expm(gt*E)
    B=[None]+[slot(t) for t in range(1,L+1)]
    res=[]
    for n in range(17,25):
        P=np.eye(N,dtype=complex)
        for k in list(range(n+1,L+1))+list(range(1,n+1)): P=B[k]@P
        G=np.linalg.inv(np.eye(N)+P); res.append(round(G[pin,pin].real,8))
    return res
print("ALF  nosymm raw1 17..24: 0.39106467 0.35351219 0.29821302 0.29821302 0.35351219 0.39106467 0.41692079 0.43497833")
print("ALF  symm   raw1 17..24: 0.40333546 0.37140904 0.32465684 0.29766374 0.32465684 0.37140904 0.40333546 0.42544558")
for symm in (False,True):
    for shift in (0,-1,1):
        print("toy symm",symm,"pulse shift",shift,run(shift,symm))
