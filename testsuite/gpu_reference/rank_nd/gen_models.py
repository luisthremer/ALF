#!/usr/bin/env python3
"""Generate models.py (ready parameter sets + ED numbers + validation record) from ed_cases.py and the validation summaries.
usage (pyALF venv): gen_models.py SUMMARY.json [SUMMARY.json ...]   -> writes models.py next to this file
The ED numbers are recomputed here (not copied from the summaries); the summaries supply the CPU-vs-ED record."""
import json
import os
import pprint
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ed_cases import CASES  # noqa: E402
from ed_validate import wfit  # noqa: E402

# GPU-test settings: dtau at which the Trotter bias is below the statistical error of the validation runs, sweeps per bin and bins
TEST = {"tV_chain4": (0.05, 400, 200), "tV_2x2_dense": (0.05, 400, 200), "tV_2x2_cb": (0.05, 400, 200),
        "Kondo_2x1_cb": (0.05, 400, 200), "Kondo_2x2_cb": (0.05, 400, 200), "SP_2x1": (0.1, 1000, 200)}
NOTES = {
    "tV_chain4": "4-site periodic chain, no duplicated bonds: Kin = <H_T>, Ener = <H> exactly.",
    "tV_2x2_dense": "2x2: every nearest-neighbour pair is a DUPLICATED bond. Checkerboard=False: propagation hops with t (Op_T '=' overwrite) but "
                    "Kin counts the pair twice (Kin = 2<H_T,prop>); Op_V is doubled. Ener is then not <H>. ED mirrors ALF.",
    "tV_2x2_cb": "2x2 with Checkerboard=True: duplicated pairs hop with 2t consistently, Kin = <H_T>, Ener = <H>.",
    "Kondo_2x1_cb": "Smallest Bilayer_square (L1=2, L2=1): duplicated c-hopping (2t, checkerboard). Ener excludes the Uf term; Pot uses Obser's JK "
                    "estimator which omits the pair terms (so Pot != <Uc term>+<H_JK>); the JK constant is JK*((2*Ns)/8 integer division) = 0.",
    "Kondo_2x2_cb": "Bilayer 2x2 (8 orbitals, 65536 states, sectors <= 4900): JK constant = JK. Same estimator remarks as Kondo_2x1_cb.",
    "SP_2x1": "Smallest Spin_Peierls lattice (Square L1=2, L2=1: 2 sites, 2 bonds = 2 phonons), SU(2) symmetric. ALF writes NO energy; "
              "Part is 2<N> (Obser counts GRC twice). Phonon bins are autocorrelated at small dtau.",
}


def main(paths):
    models = {}
    for p in paths:
        res = json.load(open(p))
        name = res["case"]
        case = CASES[name]
        sd = dict(case["sim_dict"])
        dt, nsw, nbin = TEST[name]
        ed = case["ed"](sd)
        dts = np.array(res["dtaus"])
        val = {"summary_file": os.path.basename(p), "dtaus": res["dtaus"], "nsweep": res["nsweep"], "nbin": res["nbin"], "nskip": res["nskip"],
               "run_dirs": [f"~/alf-gpu/work/t401ed/runs/<tag>/{name}_dt{d:g}" for d in res["dtaus"]], "scalars": {}, "lattice_max_abs_z": {}}
        for k in [k for k in ed if k.endswith("_scal0")]:
            ys = np.array([res["mc"][str(d)][k] for d in dts])
            es = np.array([res["mc"][str(d)][k + "_err"] for d in dts])
            if np.all(es < 1e-12):
                val["scalars"][k] = {"exact_conserved": True, "mc": float(ys[-1]), "ed": float(ed[k])}
                continue
            y0, e0, chi2 = wfit(dts, ys, es)
            i = int(np.argmin(abs(dts - dt))) if dt in dts else None
            val["scalars"][k] = {"mc_extrap": y0, "err_extrap": e0, "z_extrap": (y0 - float(ed[k])) / e0, "chi2_fit": chi2,
                                 "per_dtau": {f"{d:g}": {"mc": float(y), "err": float(e), "z": float((y - ed[k]) / e)} for d, y, e in zip(dts, ys, es)}}
            if i is not None:
                val["scalars"][k]["trotter_bias_at_test_dtau"] = float(ys[i] - ed[k])
                val["scalars"][k]["err_at_test_dtau"] = float(es[i])
        for k in [k for k in ed if k.endswith("_eqR")]:
            edv = np.asarray(ed[k]).ravel()
            vals = np.array([np.asarray(res["mc"][str(d)][k]).ravel() for d in dts])
            errs = np.array([np.asarray(res["mc"][str(d)][k + "_err"]).ravel() for d in dts])
            zl = [(wfit(dts, vals[:, j], errs[:, j])[0] - edv[j]) / wfit(dts, vals[:, j], errs[:, j])[1]
                  for j in range(edv.size) if np.all(errs[:, j] > 1e-12)]
            val["lattice_max_abs_z"][k] = float(max(abs(np.array(zl)))) if zl else None
        params = {**sd, "dtau": dt, "NSweep": nsw, "NBin": nbin}
        e_name = "Ener_scal0" if "Ener_scal0" in ed else None
        models[name] = {
            "ham_name": case["ham_name"],
            "lattice": {k: sd[k] for k in ("Lattice_type", "L1", "L2")},
            "parameters": params,
            "ed_energy": float(ed[e_name]) if e_name else None,
            "ed_energy_name": e_name,
            "ed_observables": {k: (np.asarray(v).tolist() if not np.isscalar(v) else float(v)) for k, v in ed.items()},
            "ed_eqR_layout": "array[a][b][iR]: a,b = orbital in the unit cell, iR in the order of ALF Lattice%list (pyALF Lattice.listr); "
                             "Den_eqR and Phi_eqR/Dimer_eqR are connected (Latt0 products subtracted), Green_eqR and SpinZ_eqR are not",
            "note": NOTES[name],
            "validation_cpu_vs_ed": val,
        }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models.py")
    with open(out, "w") as f:
        f.write('"""Ready parameter sets and ED references for the rank>1 / non-diagonal GPU tests (ticket 401).\n'
                'Generated by gen_models.py from ed_cases.py + the CPU validation summaries; do not edit by hand.\n'
                '`parameters` are pyALF sim_dict keys (Simulation(alf_src, ham_name, parameters, mpi=False)); dtau, NSweep, NBin are the GPU-test values.\n'
                'ED names follow the pyALF analysis output (Ener_scal0, ..., Den_eqR, ...); see ed_eqR_layout.\n"""\n\n')
        f.write("MODELS = " + pprint.pformat(models, width=130, sort_dicts=False) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:])
