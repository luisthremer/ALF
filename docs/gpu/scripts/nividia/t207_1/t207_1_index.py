#!/usr/bin/env python3
"""Ticket 207-1 (development only): index of the run directories of the T0 runs (t207_1_run.sh / t207_1_irun.sh write no provenance.json), generated AFTER
the runs from the chain logs (rc, wall, binary prefix, env) plus the binary sha256, the source state and the sha256 of run.log / info / data.h5 of every run
directory.  usage: t207_1_index.py OUT.json LOG [LOG...]   (run lines: 'TAG BIN CONFIG rc=N wall=Ts bin PREFIX env [..]')"""
import hashlib, json, os, re, socket, sys
W = os.path.expanduser('~/alf-gpu/work/t207-1')
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else None
def rd(p): return open(p).read().strip() if os.path.exists(p) else None
rows = []
for lg in sys.argv[2:]:
    for l in open(lg):
        m = re.match(r'^(?:== \S+: )?(\S+) (\S+) (\S+) rc=(\d+) wall=([0-9.]+)s bin ([0-9a-f]{16}) env \[(.*)\]', l)
        if not m: continue
        tag, b, c, rc, wall, bp, env = m.groups(); d = f'{W}/runs/{tag}/{c}'
        bs = rd(f'{W}/bin/ALF_{b}.out.sha256'); bs = bs.split()[0] if bs else None
        if bs and not bs.startswith(bp):   # a binary that was moved to bin/prefix/ is not the one that ran
            bs = None
        rows.append(dict(tag=tag, build=b, config=c, rc=int(rc), wall_s=float(wall), binary_prefix=bp, binary_sha256=bs, env=env, host=socket.gethostname(),
                         run_log_sha256=sha(d + '/run.log'), info_sha256=sha(d + '/info'), data_h5_sha256=sha(d + '/data.h5'), seeds_sha256=sha(d + '/seeds'),
                         parameters_sha256=sha(d + '/parameters'), run_dir=d, log=os.path.basename(lg)))
json.dump(dict(note='generated after the runs (t207_1_index.py); the run scripts of this ticket wrote no provenance.json; source state in srcstate',
               srcstate=rd(W + '/src/.srcstate'), runs=rows), open(sys.argv[1], 'w'), indent=1)
print(len(rows), 'runs indexed;', sum(1 for r in rows if r['binary_sha256'] is None), 'without a matching binary sha256 file')
