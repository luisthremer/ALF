#!/usr/bin/env python3
"""Ticket 207-3 (development only): write a gate spec (format of t207_d5_gate.py at gpu-port-207-2 305488ec) for genuine stride-3 cover-1
QREF runs of t207_3_chain.sh, so that the redesigned gate can be run unmodified on real data.  Not a pre-registration: bounds are inputs.

usage: t207_3_spec.py OUT.json PATH(def|log) NMAIN [--ctl N] [--flt X:N] [--fltat X@NST:NFE:N] [--gate NAME:MEASURE:POOL:BOUND[:STATISTIC] ...]
                      [--control NAME:KIND:WHAT:TRIPS[:MUSTNOT1,MUSTNOT2]] --prog-sha SHA --src-sha SHA
  main env: DEBUG=1, QHOST=1, QVAR=8, QREF=3, FROM=296 (t207_3_chain.sh); control roles: ctl (TAU_NS=0 def / TAUSWAP=8 log, no QVAR),
  flt_<x> (FAULT_TAU=eps:<x>, no QVAR), fltat_<x>_<nst>x<nfe> (plus FAULT_TAU_AT=<nst>:<nfe>).  qs 20, cover 1, refs = frozen list,
  expect_points pTT 210, emax 280, eF 280 (the measures of the gates only).
"""
import argparse, json, sys

NUL = {'def': [0, 4, 100004, 1, 100001, 2, 1608, 4808, 5, 100005], 'log': [0, 4, 100004, 1, 100001, 2, 5, 100005]}
REF = {'def': 'cpu,nvh,n4,n100004,n1,n100001,n2,n1608,n4808,n5,n100005'.split(','), 'log': 'cpu,nvh,n4,n100004,n1,n100001,n2,n5,n100005'.split(',')}
PTS = {'pTT': 210, 'emax': 280, 'eF': 280}


def mlist(nul, var):
    return ['dev'] + (['var8'] if var else []) + ['nvh', 'cpu'] + ['n%d' % c for c in nul] + ['u1', 'u2']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('path', choices=('def', 'log')); ap.add_argument('nmain', type=int)
    ap.add_argument('--ctl', type=int, default=0); ap.add_argument('--flt', action='append', default=[]); ap.add_argument('--fltat', action='append', default=[])
    ap.add_argument('--gate', action='append', default=[]); ap.add_argument('--control', action='append', default=[])
    ap.add_argument('--prog-sha', required=True); ap.add_argument('--src-sha', required=True)
    ap.add_argument('--rundir', action='append', default=[], help='ROLE=RUNDIR: binary, parameters and seeds sha256 from provenance.json (gate fix round 4)')
    a = ap.parse_args()
    p = a.path
    members = {'main': mlist(NUL[p], True)}
    env_roles = {'main': {'ALF_GPU_DEBUG': '1', 'ALF_GPU_TAUDIAG_QHOST': '1', 'ALF_GPU_TAUDIAG_QVAR': '8', 'ALF_GPU_TAUDIAG_QREF': '3',
                          'ALF_GPU_TAUDIAG_FROM': '296'}}
    n_runs = {'main': a.nmain}
    if a.ctl:
        k, v = ('ALF_GPU_TAU_NS', '0') if p == 'def' else ('ALF_GPU_TAUSWAP', '8')
        env_roles['ctl'] = {'extra_env': {k: v}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']}
        members['ctl'] = mlist(NUL[p], False); n_runs['ctl'] = a.ctl
    for f in a.flt:
        x, n = f.split(':')
        role = 'flt_' + x.replace('.', 'p').replace('-', 'm')
        env_roles[role] = {'extra_env': {'ALF_GPU_FAULT_TAU': 'eps:' + x}, 'omit_env': ['ALF_GPU_TAUDIAG_QVAR']}
        members[role] = mlist(NUL[p], False); n_runs[role] = int(n)
    for f in a.fltat:
        x, rest = f.split('@')
        nst, nfe, n = rest.split(':')
        role = 'fltat_%s_%sx%s' % (x.replace('.', 'p').replace('-', 'm'), nst, nfe)
        env_roles[role] = {'extra_env': {'ALF_GPU_FAULT_TAU': 'eps:' + x, 'ALF_GPU_FAULT_TAU_AT': '%s:%s' % (nst, nfe)},
                           'omit_env': ['ALF_GPU_TAUDIAG_QVAR']}
        members[role] = mlist(NUL[p], False); n_runs[role] = int(n)
    gates = []
    for g in a.gate:
        t = g.split(':')
        name, meas, pool, bound = t[:4]
        gates.append(dict(name=name, measure=meas, statistic=t[4] if len(t) > 4 else 'S_ens', pool=pool, bound=float(bound)))
    controls = []
    for c in a.control:
        t = c.split(':')
        d = dict(name=t[0], kind=t[1], trips=t[3])
        d['member' if t[1] == 'variant' else 'role'] = t[2]
        if len(t) > 4 and t[4]:
            d['must_not_trip'] = t[4].split(',')
        controls.append(d)
    spec = dict(spec_version=1, path=p, host_member='nvh', members=members, refs=REF[p], qs=20, nul=NUL[p], cover=1,
                expect_points={g['measure']: PTS[g['measure']] for g in gates}, prog_sha=a.prog_sha, src_sha=a.src_sha,
                env_roles=env_roles, n_runs=n_runs, gates=gates, controls=controls)
    if a.rundir:
        seeds, bins, pars = {}, set(), set()
        for rd in a.rundir:
            role, d = rd.split('=', 1)
            pr = json.load(open(d + '/provenance.json'))
            seeds.setdefault(role, []).append(pr['seeds_sha256']); bins.add(pr['binary_sha256']); pars.add(pr['parameters_sha256'])
        if len(bins) != 1 or len(pars) != 1:
            sys.exit('runs with different binaries %s or parameters %s' % (bins, pars))
        spec.update(binary_sha256=bins.pop(), parameters_sha256=pars.pop(), seeds=seeds)
        spec.update(keep=39, hostin=0, min_reprocessed=8)       # fields of the approved gate 4118a1de (KEEP 39 keeps 8 of 280 dumps)
    with open(a.out, 'x') as f:
        json.dump(spec, f, indent=1)


if __name__ == '__main__':
    main()
