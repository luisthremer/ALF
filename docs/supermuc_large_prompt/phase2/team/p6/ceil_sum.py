# summarise out/*.out of ceil_run.sh: per config, mean per-rank TFLOP/s and aggregate per tile (sum over ranks / 8 tiles, 1 tile for r1)
import sys, re, collections
for f in sys.argv[1:]:
    ntile = 1 if '/r1.' in f else 8
    d = collections.defaultdict(list); bd = collections.defaultdict(int); us = collections.defaultdict(list)
    for l in open(f):
        m = re.match(r'CEIL rank\s+(\d+) N\s*(\d+) (\S+)\s+B\s*(\d+) reps\s*(\d+) t\s*([\d.]+) TFLOPs\s*([\d.]+) bitdiff\s*(\d+)(?: usgemm\s*([\d.]+))?', l)
        if m:
            k = (int(m.group(2)), m.group(3), int(m.group(4))); d[k].append(float(m.group(7))); bd[k] += int(m.group(8)); us[k].append(float(m.group(9) or 0))
    print('==', f)
    for k in sorted(d):
        v = d[k]; print('  N %4d %-13s B %2d : ranks %2d  per-rank %6.2f TFLOP/s  per-tile aggregate %6.2f  us/gemm %8.2f  bitdiff %d' % (k[0], k[1], k[2], len(v), sum(v)/len(v), sum(v)/ntile, sum(us[k])/len(us[k]), bd[k]))
