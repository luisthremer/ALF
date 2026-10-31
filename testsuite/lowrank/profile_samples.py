#!/usr/bin/env -S python3 -I
"""Categorise the main-thread call graph of a macOS `sample` file into time fractions (self time, innermost category frame).

Usage: profile_samples.py SAMPLE_FILE [SAMPLE_FILE ...]   prints one markdown table per file.
Categories (innermost matching frame on the stack decides; ZGERU counts as rank-1 only below Upgrade2):
flush ZGEMM (lowrank_flush, submatrix_factors), measurements, hopping/interaction wraps (hop_mod, op_mmult*, op_wrap*),
stabilisation (cgr, qdrp, udv, wrapur/wrapul, mymats), ratio/border/propose (upgrade_lowrank, lowrank_*, upgrade2 other
than ZGERU, Wrapgr self time), rank-1 updates in standard (Upgrade2 -> ZGERU), other.
The raw sample files are not committed (runs/bench/).
"""
import re
import sys
import collections


def parse(path):
    lines = open(path).read().split('\n')
    i = lines.index('Call graph:') + 1
    nodes = []
    for l in lines[i:]:
        if l.startswith('Total number in stack'):
            break
        m = re.match(r'^(\s*)((?:[+!:| ]*))(\d+) (.*)$', l)
        if not m:
            continue
        col = len(m.group(1)) + len(m.group(2))
        nodes.append((col, int(m.group(3)), m.group(4).split('  (in')[0].strip(), l))
    return nodes

CATS=[ # innermost-first on the stack; (category, regex on frame name)
 ('flush ZGEMM (Lowrank_flush / Submatrix_factors)', r'lowrank_flush|submatrix_factors'),
 ('measurements', r'obser|ham_|hamiltonian|measure|tau_m|tau_p|predefined'),
 ('hopping / interaction wraps (Hop_mod, Op_mmult*, Op_wrap*)', r'hop_mod|op_mmult|op_wrap|opttypes|zdslsymm'),
 ('stabilisation (QR/UDV, CGR, wrapur/wrapul)', r'cgr|qdrp|udv|wrapur|wrapul|mymats|control_prec'),
 ('ratio / border / propose (incl. Wrapgr self time)', r'upgrade_lowrank|lowrank_|upgrade2|wrapgr'),
]
def categorise(path):
    nodes=parse(path)
    mincol=min(n[0] for n in nodes)
    heads=[k for k,n in enumerate(nodes) if n[0]==mincol]
    m=nodes[heads[0]:heads[1] if len(heads)>1 else None]
    total=m[0][1]
    # self time
    res=collections.Counter(); detail=collections.Counter()
    stack=[]
    selfs=[]
    for k,(col,cnt,name,l) in enumerate(m):
        while stack and stack[-1][0]>=col: stack.pop()
        stack.append((col,k,name))
        # children sum
        j=k+1; ch=0
        while j<len(m) and m[j][0]>col:
            if m[j][0]==col+2 or True:
                pass
            j+=1
        selfs.append(None)
    # compute self via child list
    parent=[-1]*len(m); st=[]
    for k,(col,cnt,name,l) in enumerate(m):
        while st and m[st[-1]][0]>=col: st.pop()
        parent[k]=st[-1] if st else -1
        st.append(k)
    csum=[0]*len(m)
    for k in range(len(m)):
        if parent[k]>=0: csum[parent[k]]+=m[k][1]
    for k,(col,cnt,name,l) in enumerate(m):
        s=cnt-csum[k]
        if s<=0: continue
        path=[]; p=k
        while p>=0: path.append(m[p][2].lower()); p=parent[p]
        cat='other'
        # rank-1: ZGERU under upgrade2
        if any('upgrade2' in x for x in path) and path[0].startswith('zgeru') or (any('upgrade2' in x for x in path) and any(x.startswith('zgeru') or x.startswith('cblas_zgeru') for x in path)):
            cat='rank-1 updates in standard (Upgrade2 -> ZGERU)'
        else:
            for x in path:  # innermost first
                for c,pat in CATS:
                    if re.search(pat,x): cat=c;break
                if cat!='other': break
        res[cat]+=s
        detail[(cat,m[k][2][:40])]+=s
    return total,res,detail

def main(argv):
    for p in argv[1:]:
        t, r, d = categorise(p)
        print("%s: %d samples (main thread)" % (p, t))
        print("| category | samples | fraction |\n|---|---:|---:|")
        for c, v in sorted(r.items(), key=lambda x: -x[1]):
            print("| %s | %d | %.1f %% |" % (c, v, 100.0 * v / t))
        print()


if __name__ == '__main__':
    main(sys.argv)
