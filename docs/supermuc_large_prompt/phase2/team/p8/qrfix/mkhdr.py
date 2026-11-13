#!/usr/bin/env python3
"""Generate the TEST-ONLY header copies (never the kernel source): hdr_<variant>/alf_la_sycl_qr2.hpp
  old : 7e8b6032^ (ex parity per step index s)            new : HEAD (ex parity per executed step)
  tau : new + single tau slot in the panel (plant race_tau)   v : new + single v buffer in the panel (plant race_v)
  ex  : new + single exchange slot in the apply kernel (plant race_ex)
Every copy carries the delay hooks under #ifdef QR2_TESTDELAY (without the macro the file is the original kernel text, plus plant edits):
  apply2_t: sub-group sid==1 spins QR2_DELAY iterations of a dependent FMA chain right after the exchange barrier, before reading ex
  panel2_t: sub-group of column jb-1 spins QR2_PDELAY iterations right after the step barrier, before reading tau/v (G==1: no barrier follows in the step)
The counts are read from the environment on the host at every launch (no recompile to change them)."""
import subprocess, sys, os
repo = sys.argv[1]
def git(rev): return subprocess.check_output(["git", "-C", repo, "show", rev + ":Libraries/Modules/alf_la_sycl_qr2.hpp"]).decode()
def sub(s, a, b, cnt=1):
    assert s.count(a) == cnt, (a, s.count(a)); return s.replace(a, b)
def hook(s):
    s = sub(s, "namespace alfdevqr {\n", """namespace alfdevqr {
#ifdef QR2_TESTDELAY
inline int qr2_env_int(const char *n) { const char *e = std::getenv(n); return e ? std::atoi(e) : 0; }
#endif
""")
    # apply2_t: host-side read + kernel spin
    s = sub(s, "  const int nwg = (ncg + CPW - 1) / CPW;\n  q.submit(", "  const int nwg = (ncg + CPW - 1) / CPW;\n#ifdef QR2_TESTDELAY\n  const int dly_ = G > 1 ? qr2_env_int(\"QR2_DELAY\") : 0;\n#endif\n  q.submit(")
    s = sub(s, """            it.barrier(sycl::access::fence_space::local_space);
#pragma unroll
            for (int v = 0; v < NV; ++v) {
              dot[v] = {0.0, 0.0};""", """            it.barrier(sycl::access::fence_space::local_space);
#ifdef QR2_TESTDELAY
            if (sid == 1) { double x_ = 1.0 + 1e-3 * (double)lane; for (int w_ = 0; w_ < dly_; ++w_) x_ = x_ * 1.0000001 + 1e-9; if (x_ < 0.0) ex[0] = x_; }
#endif
#pragma unroll
            for (int v = 0; v < NV; ++v) {
              dot[v] = {0.0, 0.0};""")
    # panel2_t
    s = sub(s, "  constexpr bool KEEPV = KM <= 20;\n  q.submit(", "  constexpr bool KEEPV = KM <= 20;\n#ifdef QR2_TESTDELAY\n  const int pdly_ = qr2_env_int(\"QR2_PDELAY\");\n#endif\n  q.submit(")
    s = sub(s, "        Z w = {0.0, 0.0};\n        double tr = 0.0, ti = 0.0;\n", """#ifdef QR2_TESTDELAY
        if (c == jb - 1) { double x_ = 1.0 + 1e-3 * (double)lane; for (int w_ = 0; w_ < pdly_; ++w_) x_ = x_ * 1.0000001 + 1e-9; if (x_ < 0.0) sx[0] = x_; }
#endif
        Z w = {0.0, 0.0};
        double tr = 0.0, ti = 0.0;
""")
    return s
old = hook(git("7e8b6032^")); new = hook(git("HEAD"))
v = {"old": old, "new": new}
v["tau"] = new.replace("tb[2 * buf", "tb[2 * 0"); assert v["tau"].count("tb[2 * 0") == 4
v["v"] = sub(new, "const int vo = buf * MP;", "const int vo = 0;")
v["ex"] = sub(new, "((nex & 1) * CPW * G", "(0 * CPW * G", 2)
here = os.path.dirname(os.path.abspath(__file__))
for k, s in v.items():
    d = os.path.join(here, "hdr_" + k); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "alf_la_sycl_qr2.hpp"), "w").write(s)
    print(k, len(s.splitlines()), "lines")
