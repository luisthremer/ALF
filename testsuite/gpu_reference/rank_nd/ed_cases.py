"""Cases: ALF sim_dict + ED reference with ALF output names.  Used by ed_validate.py and gen_models.py."""
import numpy as np
from edcommon import to_alf_arrays


def tV_ed(sd):
    import ed_tV
    L1, L2 = sd["L1"], sd["L2"]
    r = ed_tV.compute(L1, L2, t=sd["ham_T"], mu=sd["Ham_chem"], V=sd["ham_V"], beta=sd["beta"], checkerboard=sd["Checkerboard"])
    out = {f"{k}_scal0": r[k] for k in ("Kin", "Pot", "Part", "Ener")}
    arrs, listr = to_alf_arrays(r, L1, L2, 1)
    out.update(arrs)
    out["Den_eq0"] = np.array([r["Den_eq0[0]"]])
    return out


COMMON = {"Symm": True, "Checkerboard": False, "Ltau": 0, "Nwrap": 4}

CASES = {}
CASES["tV_chain4"] = {
    "ham_name": "tV", "ed": tV_ed,
    "sim_dict": {**COMMON, "Lattice_type": "Square", "L1": 4, "L2": 1, "ham_T": 1.0, "Ham_chem": 0.0, "ham_V": 2.0, "beta": 2.0},
}
_t22 = {**COMMON, "Nwrap": 2, "Lattice_type": "Square", "L1": 2, "L2": 2, "ham_T": 1.0, "Ham_chem": 0.0, "ham_V": 1.0, "beta": 2.0}
CASES["tV_2x2_dense"] = {"ham_name": "tV", "ed": tV_ed, "sim_dict": {**_t22, "Checkerboard": False}}
CASES["tV_2x2_cb"] = {"ham_name": "tV", "ed": tV_ed, "sim_dict": {**_t22, "Checkerboard": True}}


def Kondo_ed(sd):
    import ed_Kondo
    L1, L2 = sd["L1"], sd["L2"]
    r = ed_Kondo.compute(L1, L2, t=sd["ham_T"], mu=sd["Ham_chem"], Uc=sd["Ham_Uc"], Uf=sd["ham_Uf"], JK=sd["Ham_JK"], beta=sd["beta"], checkerboard=sd["Checkerboard"])
    out = {f"{k}_scal0": r[k] for k in ("Kin", "Pot", "Part", "Ener", "Constraint")}
    arrs, listr = to_alf_arrays(r, L1, L2, 2)
    out.update(arrs)
    out["Den_eq0"] = np.array([r["Den_eq0[0]"], r["Den_eq0[1]"]])
    return out


_K = {**COMMON, "Lattice_type": "Bilayer_square", "ham_T": 1.0, "Ham_chem": 0.0, "Ham_Uc": 1.0, "ham_Uf": 2.0, "Ham_JK": 2.0, "beta": 2.0}
CASES["Kondo_2x1_cb"] = {"ham_name": "Kondo", "ed": Kondo_ed, "sim_dict": {**_K, "L1": 2, "L2": 1, "Checkerboard": True}}
CASES["Kondo_2x1_dense"] = {"ham_name": "Kondo", "ed": Kondo_ed, "sim_dict": {**_K, "L1": 2, "L2": 1, "Checkerboard": False}}
CASES["Kondo_2x2_cb"] = {"ham_name": "Kondo", "ed": Kondo_ed, "sim_dict": {**_K, "L1": 2, "L2": 2, "Checkerboard": True}}


def SP_ed(sd):
    import ed_SpinPeierls as E
    L1, L2 = sd["L1"], sd["L2"]
    r = E.compute(L1, L2, J=sd["Ham_Jx"], U=sd["Ham_U"], lam=sd["Ham_Lambda"], om0=sd["Ham_Omega0"], beta=sd["beta"], nmax=sd.get("_nmax", 12))
    out = {"Part_scal0": 2.0 * r["Part"],  # Obser: ZRho = N_SUN*(sum GRC + sum GRC): counts every electron twice
            "Pot_scal0": r["Pot"], "PhiXY_scal0": r["PhiXY[0]"], "PhiXY_scal1": r["PhiXY[1]"]}
    arrs, listr = E.alf_arrays(r, L1, L2)
    out.update(arrs)
    return out


CASES["SP_2x1"] = {"ham_name": "Spin_Peierls", "ed": SP_ed,
                   "sim_dict": {"Symm": False, "Checkerboard": False, "Ltau": 0, "Nwrap": 4, "Lattice_type": "Square", "L1": 2, "L2": 1,
                                "Ham_Jx": 1.0, "Ham_Jy": 0.0, "Ham_U": 4.0, "Ham_h": 0.0, "Ham_g_factor": 2.0, "Ham_Lambda": 0.1,
                                "Ham_Omega0": 2.0, "beta": 2.0}}
