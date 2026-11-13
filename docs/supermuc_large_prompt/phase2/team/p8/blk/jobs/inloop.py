import sys,os
for d in sys.argv[1:]:
    f=os.path.join(d,'run.log')
    if not os.path.exists(f): continue
    rows=[]
    for l in open(f):
        if 'TEAM_TIMERS' in l and '[s]:' in l:
            x=l.split('[s]:')[-1].split()
            if len(x)>=8: rows.append([float(v) for v in x])
    if not rows: continue
    n=rows[0][0]
    m=lambda i: sum(r[i] for r in rows)/len(rows)/n
    print("%-24s ranks %2d sweep %.3f upd %.3f stab %.3f cgr %.3f obs %.3f"%(os.path.basename(d),len(rows),m(1),m(2)+m(3),m(4)+m(5),m(6),m(7)))
