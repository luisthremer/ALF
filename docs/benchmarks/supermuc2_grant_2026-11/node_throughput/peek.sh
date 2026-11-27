#!/bin/bash
# print per running run: name, number of bursts so far (bin writes) and seconds since start / since last write
now=$(date +%s%N)
for d in /dev/shm/nt_*/cpu_*; do
  [ -f $d/watch.txt ] || continue
  python3 - "$d" "$now" <<'PY'
import sys,os
d,now=sys.argv[1],int(sys.argv[2])
ev=sorted({int(l.split()[0]) for l in open(d+"/watch.txt") if len(l.split())==3})
st=os.stat(d+"/parameters").st_mtime_ns
bur=[];last=None
for t in ev:
    if last is None or t-last>2e9: bur.append([t,t])
    else: bur[-1][1]=t
    last=t
ends=[b[1] for b in bur[1:]]
gaps=[(b-a)/1e9 for a,b in zip(ends,ends[1:])]
print(os.path.basename(d), "bins:",len(ends), "since start: %.0f s"%((now-st)/1e9),
      "first bin at: %s"%("%.0f s"%((ends[0]-st)/1e9) if ends else "-"), "sweep gaps:",["%.1f"%g for g in gaps])
PY
done
