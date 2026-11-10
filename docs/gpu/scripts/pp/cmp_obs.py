import h5py, numpy as np, sys, os
R=os.path.expanduser('~/alf-gpu/work/pp/runs/')
def load(tag,i):
    f=h5py.File(R+f'bench_{tag}_n3_p25_d64/Temp_{i}/data.h5','r'); return f
f=load('base',0); 
names=[]; f.visititems(lambda n,o: names.append(n) if isinstance(o,h5py.Dataset) else None)
print([n for n in names if 'SpinZ_eq' in n or 'Ener' in n or 'Den_eq' in n][:12])
for i,U in enumerate((3.0,5.25,8.0)):
    for tag in ('opt','optnw20','optall','base2'):
        try: g=load(tag,i)
        except Exception as e: print(tag,'missing'); continue
        out=[]
        for key in ('SpinZ_eq/obser','Ener_scal/obser','Den_eq/obser','Green_eq/obser'):
            try:
                a=np.array(load('base',i)[key]); b=np.array(g[key])
                out.append(f'{key.split("/")[0]}: max|d|={np.max(np.abs(a-b)):.2e} (max|a|={np.max(np.abs(a)):.2e})')
            except Exception as e: out.append(f'{key.split("/")[0]}: n/a')
        print(f'U={U} {tag}:', ' | '.join(out))
