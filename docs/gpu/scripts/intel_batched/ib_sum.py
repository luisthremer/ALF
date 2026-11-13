#!/usr/bin/env python3
"""Summarise one ib_s1 output directory: python3 ib_sum.py <out/JOBID>  -> markdown tables on stdout."""
import collections, glob, os, re, sys

d = sys.argv[1]


def kv(line):
    return dict(re.findall(r'(\w+)=(\S+)', line))


def lines(f, key):
    if not os.path.exists(f):
        return []
    return [kv(l) for l in open(f) if l.startswith(key + ' ')]


print('## triad (one process, tile alone)\n\n| n | MB/array | us | GB/s |\n|--|--|--|--|')
for r in lines(f'{d}/r1.txt', 'triad'):
    print(f"| {r['n']} | {r['MB_per_array']} | {r['t_us']} | {r['GBs']} |")

print('\n## launch latency (SYCL, in-order queue)\n\n| file | sync us | async submit us | async total us | copy sync us | copy async us |\n|--|--|--|--|--|--|')
for f in sorted(glob.glob(f'{d}/launch*.txt')) + [f'{d}/r1.txt']:
    for r in lines(f, 'launch'):
        print(f"| {os.path.basename(f)} | {r['sync_us']} | {r['async_submit_us']} | {r['async_total_us']} | {r['copy_sync_us']} | {r['copy_async_us']} |")
for f in sorted(glob.glob(f'{d}/omp*.txt')):
    print(f'\n{os.path.basename(f)}: ' + open(f).read().strip())


def gemmtab(files, title):
    acc = collections.defaultdict(list)
    for f in files:
        for r in lines(f, 'gemm'):
            acc[(int(r['N']), int(r['K']), int(r['B']), r['mode'])].append(float(r['TFs']))
    if not acc:
        return
    print(f'\n## gemm: {title} (TFLOP/s per process; sum over processes in brackets)\n')
    print('| N | K | B | loop | batch |\n|--|--|--|--|--|')
    keys = sorted({k[:3] for k in acc})
    for k in keys:
        cells = []
        for m in ('loop', 'batch'):
            v = acc.get(k + (m,))
            cells.append('' if not v else f'{sum(v) / len(v):.2f} [{sum(v):.2f}]')
        print(f'| {k[0]} | {k[1]} | {k[2]} | {cells[0]} | {cells[1]} |')


gemmtab([f'{d}/r1.txt'], '1 process, tile alone')
for P in (2, 4):
    for c in (1, 4):
        gemmtab(sorted(glob.glob(f'{d}/p{P}_ccs{c}_*.txt')), f'{P} processes on tile 0, CCS {c}')
gemmtab(sorted(glob.glob(f'{d}/n8_t*.txt')), '8 tiles busy, 1 process each')


def blktab(files, title):
    acc = collections.defaultdict(list)
    for f in files:
        for r in lines(f, 'block'):
            acc[(int(r['N']), int(r['W']), r['mode'])].append(float(r['us_per_walker_block']))
    if not acc:
        return
    print(f'\n## update-block emulation: {title} (us per walker-block, mean over processes)\n')
    print('| N | W | today (W x single, wait each) | batched (1 launch/op, 1 wait) | batched, no wait |\n|--|--|--|--|--|')
    for k in sorted({k[:2] for k in acc}):
        c = [acc.get(k + (m,)) for m in ('today', 'batched', 'batched_nowait')]
        print(f'| {k[0]} | {k[1]} | ' + ' | '.join('' if not v else f'{sum(v) / len(v):.1f}' for v in c) + ' |')


blktab([f'{d}/r1.txt'], '1 process, tile alone')
blktab(sorted(glob.glob(f'{d}/blk4_ccs4_*.txt')), '4 processes per tile, CCS 4 (production layout)')
for f in ('wtest.txt', 'wtest_plant.txt'):
    if os.path.exists(f'{d}/{f}'):
        print(f'\n{f}:\n```\n' + open(f'{d}/{f}').read().strip() + '\n```')
