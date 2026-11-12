import sqlite3, glob, sys
iv=[]; spans=[]
for f in sorted(glob.glob('rank*.sqlite')):
    c=sqlite3.connect(f); r=c.execute("select min(start),max(end),count(*) from CUPTI_ACTIVITY_KIND_KERNEL").fetchone(); spans.append((f,)+r)
    iv+= c.execute("select start,end from CUPTI_ACTIVITY_KIND_KERNEL").fetchall()
for s in spans: print(s[0], round((s[1]-spans[0][1])/1e9,2), round((s[2]-spans[0][1])/1e9,2), s[3])
iv.sort(); busy=0; cs,ce=iv[0]
for s,e in iv[1:]:
    if s>ce: busy+=ce-cs; cs,ce=s,e
    else: ce=max(ce,e)
busy+=ce-cs
t0=min(s[1] for s in spans); t1=max(s[2] for s in spans)
tot=sum(e-s for s,e in iv)
print("window s", (t1-t0)/1e9, "GPU busy union s", busy/1e9, "frac", busy/(t1-t0), "sum of kernel durations s", tot/1e9, "mean concurrency while busy", tot/busy)
