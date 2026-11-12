import csv, sys, collections, glob
root = sys.argv[1]
D = {}   # dispatch -> dict
for p in sorted(glob.glob(root + "/pass*/prof/r_counter_collection.csv")):
    for r in csv.DictReader(open(p)):
        k = int(r["Dispatch_Id"]); d = D.setdefault(k, {"name": r["Kernel_Name"], "t": int(r["End_Timestamp"]) - int(r["Start_Timestamp"]),
            "grid": int(r["Grid_Size"]), "wg": int(r["Workgroup_Size"]), "lds": int(r["LDS_Block_Size"]), "scr": int(r["Scratch_Size"]),
            "vgpr": int(r["VGPR_Count"]), "c": {}})
        d["c"][r["Counter_Name"]] = d["c"].get(r["Counter_Name"], 0) + float(r["Counter_Value"])
def grp(n):
    l = n.lower()
    if "tensile" in l or "cijk" in l: return "rocBLAS gemm (Tensile)"
    if "rocsolver" in l or "getrf" in l or "laswp" in l or "trsm" in l: return "rocSOLVER/trsm"
    if "k_qr_panel" in l: return "QR panel"
    if "k_apply_left" in l: return "QR apply"
    if "omp" in l or "_omp_fn" in l or "main._" in l or "gpu_" in l: return "OpenMP-offload kernels"
    return "HIP kernels (" + n[:30] + ")"
G = collections.defaultdict(lambda: collections.defaultdict(float))
N = collections.defaultdict(lambda: collections.defaultdict(float))
for k, d in D.items():
    g = grp(d["name"]); t = d["t"]; G[g]["n"] += 1; G[g]["t"] += t
    for key in ("OccupancyPercent", "MeanOccupancyPerCU", "VALUBusy", "MemUnitBusy", "L2CacheHit", "MemUnitStalled"):
        if key in d["c"]: G[g][key] += d["c"][key] * t
    if "FETCH_SIZE" in d["c"]: G[g]["bytes"] += (d["c"]["FETCH_SIZE"] + d["c"].get("WRITE_SIZE", 0)) * 1024
    for key in ("wg", "grid", "vgpr", "lds", "scr"): G[g][key] += d[key] * t
tot = sum(g["t"] for g in G.values())
print("total kernel time %.2f s, %d dispatches" % (tot/1e9, len(D)))
print("%-34s %7s %8s %6s | %5s %5s %5s %5s %5s | %7s %5s | %4s %4s %5s %5s" % ("group", "n", "ms", "%", "occ%", "VALU", "MemB", "L2hit", "stall", "GB/s", "%peak", "wg", "vgpr", "lds", "scr"))
for g, v in sorted(G.items(), key=lambda x: -x[1]["t"]):
    t = v["t"]
    if t / tot < 0.004: continue
    bw = v["bytes"] / t   # B/ns = GB/s
    print("%-34s %7d %8.0f %6.1f | %5.1f %5.1f %5.1f %5.1f %5.1f | %7.0f %5.1f | %4.0f %4.0f %5.0f %5.0f" % (g[:34], v["n"], t/1e6, 100*t/tot,
        v["OccupancyPercent"]/t, v["VALUBusy"]/t, v["MemUnitBusy"]/t, v["L2CacheHit"]/t, v["MemUnitStalled"]/t, bw, 100*bw/1638.4,
        v["wg"]/t, v["vgpr"]/t, v["lds"]/t, v["scr"]/t))
