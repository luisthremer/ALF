"""All serial scenarios of test_branch_parameters.json with the GPU=host build in check mode: GPU line + CGR check per scenario."""
import json, os, sys, re
build = sys.argv[1]
root = os.path.expanduser("~/work/ALF_P2/t5/src_" + build)
exe = os.path.join(root, "Prog/ALF.out")
sys.path.insert(0, os.path.join(root, "testsuite/Lowrank.tests"))
from compare_schemes import make_text, prepare, execute
start = os.path.join(root, "Scripts_and_Parameters_files/Start")
template = open(os.path.join(start, "parameters")).read()
d = json.load(open(os.path.join(root, "testsuite/test_branch_parameters.json")))
os.environ.update(ALF_GPU_DEBUG="1", ALF_GPU_CGR_LU="2", OMP_NUM_THREADS="1")
for name, (model, pars) in d.items():
    if isinstance(pars, list): print("%-6s scen-%-22s SKIP tempering list" % (build, name)); continue
    text, missing = make_text(template, pars, 2, None)
    if missing: print("%-6s scen-%-22s SKIP missing %s" % (build, name, missing)); continue
    w = os.path.expanduser("~/work/ALF_P2/t5/runs/%s/scen-%s" % (build, name))
    prepare(w, [text], os.path.join(start, "seeds"), 1, "serial")
    rc = execute(exe, w, None, 1)
    txt = open(os.path.join(w, "info"), errors="replace").read() if os.path.exists(os.path.join(w, "info")) else ""
    act = re.search(r"GPU path \(gpu_state\): (.*)", txt)
    chk = re.findall(r"CGR calls compared (\d+), max \|phase\(tracked\) - phase\(LU\)\| +(\S+)", txt)
    print("%-6s scen-%-22s rc=%s %s | check %s" % (build, name, rc, act.group(1)[:90] if act else "no GPU line", chk or "none"))
