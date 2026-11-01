#!/usr/bin/env python3
"""Ticket 207-2: compare the new pipeline's qref.txt of the kept dev5 dumps (t207_d5_repro.sh) with the dev5 record of the same dumps.

usage: t207_d5_repro.py NEWWORK TAG [TAG...]     (reads NEWWORK/runs/TAG/hub12x12_cb_p/qref.txt and dev5/runs/TAG/hub12x12_cb_p/qref.txt)
The dev5 record is in the old format (ES10.3, 5 numbers per method for dump version 1, numbers of three-digit exponents without 'E'); it is
read here with a deliberately lenient reader that is NOT part of the pipeline.  The new file is read by the strict parser.
Per run: points compared (the kept dumps, by dump number), members, numbers; numbers that agree after rounding the new value to ES10.3
('equal'), that differ by one unit in the last printed digit of the old format ('1 ulp'), or by more ('BAD'); header fields that differ.
A number counts as reproduced if it is 'equal' or '1 ulp' where the new value's last printed digit is a 5 (the new value has 5 digits, the old 4:
rounding the new value again to 4 digits can differ from the old direct rounding only then).  Exit status 1 if any BAD number, header or member difference exists.
"""
import importlib.util, math, os, re, sys
sys.dont_write_bytecode = True      # the loaders below import sibling scripts by path: no __pycache__ in the scripts directory

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('qp', os.path.join(HERE, 't207_d5_qparse.py'))
QP = importlib.util.module_from_spec(spec)
spec.loader.exec_module(QP)
OLD = '/home/' + os.environ.get('USER', 'luis') + '/alf-gpu/work/t207/dev5'
TOK = re.compile(r'[-+]?\d\.\d{3}E[-+]\d{2,3}|\d\.\d{3}[-+]\d{3}')


def read_old(fn):
    pts = {}
    for line in open(fn):
        if not line.startswith('QREF'):
            continue
        segs = line.rstrip('\n').split(' | ')
        h = segs[0].split()
        hd = {h[i]: h[i + 1] for i in range(1, len(h) - 1, 2)}
        m = {}
        names = []
        for s in segs[1:]:
            t = s.split()
            m[t[0]] = TOK.findall(s[len(t[0]):])
            names.append(t[0])
        pts[int(hd['file'])] = (hd, m, names)
    return pts


def main():
    work = sys.argv[1]
    bad_total = 0
    print('%-12s %6s %8s %8s %8s %8s %6s' % ('run', 'points', 'numbers', 'equal', '1 ulp', 'BAD', 'other'))
    for tag in sys.argv[2:]:
        new = QP.load('%s/runs/%s/hub12x12_cb_p/qref.txt' % (work, tag), allow_v1=tag.startswith('s_'))
        old = read_old('%s/runs/%s/hub12x12_cb_p/qref.txt' % (OLD, tag))
        nnum = neq = nulp = nbad = nother = 0
        worst = 0.0
        for p in new.points:
            k = p.hd['file']
            if k not in old:
                nother += 1; print('  %s: dump %d not in the dev5 record' % (tag, k)); continue
            hd, m, names = old[k]
            for key in ('pass', 'nst', 'nfe', 'br', 'hbr'):
                if int(hd[key]) != p.hd[key]:
                    nother += 1; print('  %s dump %d: header %s differs' % (tag, k, key))
            if names != p.names:
                nother += 1; print('  %s dump %d: member list differs %s vs %s' % (tag, k, names, p.names)); continue
            for nm in names:
                ot = m[nm]
                nv = p.m[nm]
                if len(ot) > len(nv) or (len(ot) != len(nv) and len(ot) != 5):
                    nother += 1; print('  %s dump %d member %s: %d old numbers, %d new' % (tag, k, nm, len(ot), len(nv))); continue
                for t, x in zip(ot, nv):
                    nnum += 1
                    if 'E' not in t:
                        nother += 1; print('  %s dump %d member %s: old token without E: %s' % (tag, k, nm, t)); continue
                    ov = float(t)
                    if x == ov == -1.0:
                        neq += 1; continue
                    r = float('%.3E' % x)
                    if r == ov:
                        neq += 1; continue
                    u = 10.0 ** (math.floor(math.log10(abs(ov))) - 3)
                    d = abs(r - ov) / u
                    worst = max(worst, d)
                    if abs(d - 1.0) < 0.01 and ('%.4E' % x).split('E')[0][-1] == '5':      # double rounding: the 5th digit is a 5
                        nulp += 1
                    else:
                        nbad += 1; print('  %s dump %d member %s: old %s new %.4E' % (tag, k, nm, t, x))
        print('%-12s %6d %8d %8d %8d %8d %6d   (%d of the %d new points in the old record)' % (tag, len(new.points), nnum, neq, nulp, nbad, nother,
              sum(1 for p in new.points if p.hd['file'] in old), len(new.points)))
        bad_total += nbad + nother
    sys.exit(1 if bad_total else 0)


if __name__ == '__main__':
    main()
