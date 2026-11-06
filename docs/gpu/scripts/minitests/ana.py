import sys,math
# mean and standard error of a scalar observable file (value per bin, first column = count, then (re,im) sign)
def load(path):
    v=[];sg=[]
    for l in open(path):
        l=l.replace('(',' ').replace(')',' ').replace(',',' ')
        x=l.split()
        if len(x)>=4: v.append(float(x[1])); sg.append(float(x[3]))
    return v,sg
def stat(path,skip=0):
    v,s=load(path); v=v[skip:]; n=len(v)
    m=sum(v)/n; var=sum((a-m)**2 for a in v)/(n-1); return m, math.sqrt(var/n), n
if __name__=="__main__":
    ref=float(sys.argv[2]) if len(sys.argv)>2 else None
    m,e,n=stat(sys.argv[1], int(sys.argv[3]) if len(sys.argv)>3 else 0)
    print("%s: %.5f +- %.5f (%d bins)"%(sys.argv[1].split('/')[-2],m,e,n)+("   dev from ref %.2f sigma"%(abs(m-ref)/e) if ref is not None else ""))
