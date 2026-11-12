#!/usr/bin/env python3
"""Per-call time of the cuSOLVER calls (geqrf, ormqr/unmqr, getrf) in an nsys sqlite export of an ALF GPU run (speedup S4).
usage: nsys_solver_calls.py PROFILE.sqlite
Kernels are ordered by start time and cut into runs of consecutive library kernels (a run ends at an OpenMP kernel `nvkernel_*`
or at cuSOLVER's `copy_info_kernel`; a `setup_vhat_l_kernel` after non-ormqr kernels starts a new run).  A run containing
`geqr2_smem_domino_fast` is a geqrf, one containing `setup_vhat_l_kernel` / `larft_T_64_batch` an ormqr, one containing
`getrf_pivot` a getrf (cuSOLVER's recursive LU: its own trsm/gemm kernels are in the run).  wall = first kernel start to
last kernel end of the run (GPU timeline, includes gaps), ker = sum of kernel durations.  Read-only on the sqlite file."""
import sqlite3, sys, collections
c = sqlite3.connect("file:%s?mode=ro" % sys.argv[1], uri=True)
rows = c.execute("select k.start,k.end,s.value,k.gridX from CUPTI_ACTIVITY_KIND_KERNEL k join StringIds s on s.id=k.shortName "
                 "order by k.start").fetchall()
span = (rows[-1][1] - rows[0][0]) / 1e9
QRF = {'geqr2_smem_domino_fast'}
ORM = {'setup_vhat_l_kernel', 'larft_T_64_batch'}
LU = {'getrf_pivot', 'ipiv_lower_small', 'create_pivot_v2', 'ipiv_lower_diag', 'ipiv_64_to_32_ker'}
runs = []; cur = []
for r in rows:
    n = r[2]
    if n.startswith('nvkernel') or n == 'copy_info_kernel':
        if cur: runs.append(cur); cur = []
        continue
    if cur and n in ORM and not ({x[2] for x in cur} & ORM):
        runs.append(cur); cur = []
    cur.append(r)
if cur: runs.append(cur)
agg = collections.defaultdict(lambda: [0, 0.0, 0.0, []])
for run in runs:
    names = {r[2] for r in run}
    if names & QRF: k = 'geqrf'
    elif names & ORM: k = 'ormqr'
    elif names & LU: k = 'getrf'
    else: k = 'other:' + ','.join(sorted(names))[:60]
    wall = (run[-1][1] - run[0][0]) / 1e6; ks = sum(r[1] - r[0] for r in run) / 1e6
    a = agg[k]; a[0] += 1; a[1] += wall; a[2] += ks; a[3].append(wall)
tot = sum(r[1] - r[0] for r in rows) / 1e6
print("profile span %.2f s, total kernel time %.1f ms" % (span, tot))
for k, a in sorted(agg.items(), key=lambda x: -x[1][1]):
    w = sorted(a[3]); med = w[len(w) // 2]
    print("%-45s n=%6d wall_sum=%9.1f ms ker_sum=%9.1f ms median/call=%6.2f ms max=%6.2f ms"
          % (k[:45], a[0], a[1], a[2], med, w[-1]))
