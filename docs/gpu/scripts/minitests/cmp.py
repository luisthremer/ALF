import sys
def vals(p):
    r=[]
    for l in open(p):
        l=l.replace('(',' ').replace(')',' ').replace(',',' ').split()
        if len(l)>=3: r.append(float(l[1]))
    return r
a=vals(sys.argv[1]); b=vals(sys.argv[2])
if not a or len(a)!=len(b): print("NA"); sys.exit()
m=max(abs(x-y)/max(1e-30,abs(x)) for x,y in zip(a,b)); print("%.1e"%m)
