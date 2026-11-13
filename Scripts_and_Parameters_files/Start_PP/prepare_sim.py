#!/usr/bin/env python3
"""Set up a parameter scan for ALF's parallel-parameters mode (PARALLEL_PARAMS, "PP") or for parallel tempering.

The template `parameters` file may contain, instead of a single value,
  * a set   `ham_U = {0.5d0, 1.5d0, 3.0d0}`               -> one simulation per value;
  * a linked list `L2 = [4, 6, 8; VAR_lattice, L1]`      -> varies together with VAR_lattice/L1
    (same length as the set of L1, element i goes with element i of L1).
All independent sets are combined as a Cartesian product. For every combination Temp_<i>/parameters is written
(the template with the scalar values patched in) plus a symbolic link named after the varied values
(e.g. ham_U_1.5d0_l1_6_... -> Temp_3). Groups listed in RULES are only varied when their condition holds
(e.g. sets in VAR_tV are ignored unless ham_name = "tV"); otherwise their template text is left as it is.

Afterwards: put `seeds` (one integer per MPI rank) next to `parameters`, set mpi_per_parameter_set in &VAR_TEMP
and run `mpiexec -n <number of Temp_i * mpi_per_parameter_set> $ALF_DIR/Prog/ALF.out` (ALF compiled with PP or
TEMPERING). See README.

Requires f90nml (pip install f90nml; it is also a dependency of pyALF).
"""
import argparse
import f90nml
import itertools
import os
import subprocess
from collections import defaultdict
import json
import re

#TODO: Caps sensitivity in keys has to be removed!
#TODO: Check whether second key exists...

def extract_lists(outer_key_list: list, current_key: str, some_dirty_list: list):
    outer_key_list = [str(key).lower() for key in outer_key_list]
    if type(some_dirty_list) != list:
        raise ValueError("Type list is expected.")
    
    dirty_copy = list(some_dirty_list)
    def extract_connection(dirty_copy: list):
        if dirty_copy[0]!="[":
            raise ValueError(f"Expected first element to be '[' but obtained {dirty_copy[0]}")
        if dirty_copy[-1]!="]":
            raise ValueError(f"Expected last element to be ']' but obtained {dirty_copy[-1]}")
        if str(dirty_copy[-3]).lower() not in outer_key_list:
            raise KeyError(f"First connection identifier for current key {current_key} with value {dirty_copy} has to be a namelist identifier. Choose from {outer_key_list}.")
            
        if dirty_copy[-4]!=";":
            raise ValueError(f"Expected ';' as separator but obtained {dirty_copy[-4]}")
        else:
            inner_key = dirty_copy[-2]
            outer_key = dirty_copy[-3]
            dirty_copy.pop(-4)
            dirty_copy.pop(-3)
            dirty_copy.pop(-2)
            return dirty_copy, outer_key, inner_key

    if "[" not in dirty_copy:
        return [el for el in dirty_copy if "{" !=el and "}"!=el],(None, None)
    else:
        dirty_copy, outer_key, inner_key = extract_connection(dirty_copy)
        return [el for el in dirty_copy if "[" !=el and "]"!=el], (outer_key, inner_key)

def cleanup_dirty_list(nml):
    keys1=nml.keys()
    extracted_information={}

    for key in keys1:
        if key == "var_errors":
            continue
        if key == "var_max_stoch":
            continue

        if key not in extracted_information:
            extracted_information[key] = {}

        for key2 in nml[key].keys():
            nml_el = nml[key][key2]
            current_key=(key,key2)
            if type(nml_el) == list:
                extracted_information[key][key2] = extract_lists(keys1, current_key, nml_el)
    return extracted_information

def apply_patch(original_nml_path, patch_nml, save_path):
    f90nml.patch(original_nml_path, patch_nml, save_path)

def flatten_dict(data):
    flat_data = {}
    for outer_k, inner_dict in data.items():
        for inner_k, v in inner_dict.items():
            if isinstance(v, tuple) and len(v) == 2:
                lst, link = v
                
                if isinstance(link, tuple) and len(v)==2:
                    outer_key, inner_key = link
                    if link == (None, None):
                        link = link
                    elif not isinstance(outer_key, str):
                        raise TypeError(f"outer_key has to be of type str but got {type(outer_key)}")
                    elif not isinstance(inner_key, str):
                        raise TypeError(f"inner key has to be of type str but got {type(inner_key)}.")
                    
                    if isinstance(outer_key,str) and isinstance(inner_key, str):
                        link = (outer_key.lower(), inner_key.lower())
                    if type(inner_key)!=type(outer_key):
                        raise ValueError(f"Inner and outer key are expected to be the same type, but have type (outer, inner): {type(outer_key)}{type(inner_key)} instead.")
                else:
                    raise TypeError(f"Link not successfully set.")
                
                super_key = (outer_k, inner_k)
                lst = lst if isinstance(lst, list) else [lst]
                flat_data[super_key] = (lst, link)
    return flat_data

def group_something(flat_data):
    groups = defaultdict(list)

    for key, (_, link) in flat_data.items():
        if link == (None, None):
            master = key
        else:
            if link not in flat_data:
                raise ValueError(f"Link '{link}' not found in any dictionary keys.")
            master = link
        groups[master].append(key)
    return groups

def zip_linked(flat_data, groups):
    group_combos = []
    for master, keys_in_group in groups.items():
        lengths = [len(flat_data[k][0]) for k in keys_in_group]
        temp=[flat_data[k][0] for k in keys_in_group]
        if len(set(lengths)) > 1:
            raise ValueError(f"Length mismatch in linked group: {keys_in_group}. For master_key {master}.")
        list_len = lengths[0]
        zipped_list = [{k: flat_data[k][0][i] for k in keys_in_group} for i in range(list_len)] # For each key group (flattend key), for each element i append the element. All combinations, not filtered yet.
        group_combos.append(zipped_list)
    return group_combos

def cartesian_product(group_combos):
    raw_combos = []
    for combo_tuple in itertools.product(*group_combos):
        merged = {}
        for sub_dict in combo_tuple:
            merged.update(sub_dict)
        raw_combos.append(merged)
    return raw_combos

def unflatten(data, raw_combos, conditions, defaults=None):
    final_results = []
    for flat_combo in raw_combos:
        nested = defaultdict(dict)
        for (outer_k, inner_k), val in flat_combo.items():
            nested[outer_k][inner_k] = val
            
        for outer_k, inner_dict in data.items():
            if not inner_dict:
                nested[outer_k] = {}

        for target_dict, (check_key, required_val) in conditions.items():
            actual_val = None
            for outer, inner_dict in nested.items():
                if check_key in inner_dict:
                    actual_val = inner_dict[check_key]
                    break
            
            if actual_val is None and defaults:
                # Condition parameter not varied: take its scalar value from the template.
                actual_val = defaults.get(check_key)
            if actual_val != required_val:
                if target_dict in nested:
                    del nested[target_dict]

        final_results.append(dict(nested))
    return final_results

def safely_deduplicate(final_results):
    seen = set()
    unique_results = []
    for d in final_results:
        dict_str = json.dumps(d, sort_keys=True)
        if dict_str not in seen:
            seen.add(dict_str)
            unique_results.append(d)
    return unique_results

def _sanitize_for_name(value):
    text = str(value).strip()
    text = re.sub(r"\s+", "-", text)
    text = re.sub(r"[^A-Za-z0-9_.-]", "-", text)
    text = text.strip("-_.")
    return text or "na"

def build_link_name(patch):
    parts = []
    for outer_key in sorted(patch.keys()):
        inner_dict = patch[outer_key]
        for parameter, value in sorted(inner_dict.items()):
            parts.append(f"{_sanitize_for_name(parameter)}_{_sanitize_for_name(value)}")
    return "_".join(parts) if parts else "sim"

# Groups that are kept only if a parameter has a given value (group -> (parameter, value)).
RULES = {
    'var_hubbard': ('ham_name', 'Hubbard'),
    'var_tv': ('ham_name', 'tV')
}


def main(original_param_path, rules, outdir=".", links=True):
    """Write outdir/Temp_<i>/parameters for every combination; return the list of patches (one dict per Temp_i)."""
    original_nml = f90nml.read(original_param_path)
    extracted_information = cleanup_dirty_list(original_nml)
    data = extracted_information

    flattend = flatten_dict(data)
    grouped = group_something(flatten_dict(data))
    linked = zip_linked(flattend, grouped)
    defaults = {k: v for group in original_nml.values() for k, v in group.items() if not isinstance(v, list)}
    l_of_patches = safely_deduplicate(unflatten(data, cartesian_product(linked), rules, defaults))

    for i, el in enumerate(l_of_patches):
        sim_dir = f"Temp_{i}"
        os.makedirs(os.path.join(outdir, sim_dir), exist_ok=True)
        apply_patch(original_param_path, el, os.path.join(outdir, sim_dir, "parameters"))

        if links:
            link_name = os.path.join(outdir, build_link_name(el))
            if os.path.lexists(link_name):
                os.unlink(link_name)
            os.symlink(os.path.join(".", sim_dir), link_name)
    return l_of_patches


SELF_TEST_TEMPLATE = """&VAR_ham_name
ham_name = {"Hubbard"}
/
&VAR_lattice
L1           = {4, 6, 8}
L2           = [4, 6, 8; VAR_lattice, L1]
Lattice_type = "Square"
Model        = "Hubbard"
/
&VAR_Model_Generic
N_Phi        = [8, 18, 32; VAR_lattice, L1]
Beta         = [4.d0, 6.d0, 8.d0; VAR_lattice, L1]
Dtau         = 0.1d0
/
&VAR_TEMP
mpi_per_parameter_set = 2
/
&VAR_Hubbard
ham_U        = {1.0d0, 4.0d0}
ham_T        = 1.d0
/
&VAR_tV
ham_V        = {1.d0, 2.d0}
/
&VAR_errors
n_skip = 1
/
"""


def plain_namelist_check(path, skip=()):
    """Standard-library check: balanced &group ... / blocks, one 'key = value' per line, no set/link syntax left
    (except in the groups in skip)."""
    groups, open_group = [], None
    for line in open(path):
        line = line.split("!")[0].strip()
        if not line:
            continue
        if line.startswith("&"):
            assert open_group is None, f"{path}: nested group {line}"
            open_group = line[1:].lower()
        elif line == "/":
            assert open_group is not None, f"{path}: '/' outside a group"
            groups.append(open_group)
            open_group = None
        else:
            assert open_group is not None and "=" in line, f"{path}: bad line {line!r}"
            assert open_group in skip or not any(c in line for c in "{}[];"), f"{path}: set/link syntax left: {line!r}"
    assert open_group is None, f"{path}: group {open_group} not closed"
    return groups


def self_test():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        tmpl = os.path.join(d, "parameters")
        open(tmpl, "w").write(SELF_TEST_TEMPLATE)
        patches = main(tmpl, RULES, outdir=d)
        # VAR_tV/ham_V is a set too, but ham_name = "Hubbard": it must not multiply the number of sets.
        assert len(patches) == 6, f"expected 3 lattices x 2 U = 6 sets, got {len(patches)}"
        seen = set()
        for i in range(len(patches)):
            p = os.path.join(d, f"Temp_{i}", "parameters")
            plain_namelist_check(p, skip=("var_tv",))
            nml = f90nml.read(p)
            l1, l2 = nml["var_lattice"]["l1"], nml["var_lattice"]["l2"]
            assert l1 == l2 and nml["var_model_generic"]["n_phi"] == l1 * l1 // 2, f"{p}: linked values wrong"
            assert nml["var_model_generic"]["beta"] == float(l1) and nml["var_model_generic"]["dtau"] == 0.1
            assert nml["var_errors"]["n_skip"] == 1 and nml["var_temp"]["mpi_per_parameter_set"] == 2
            seen.add((l1, nml["var_hubbard"]["ham_u"]))
            link = os.path.join(d, build_link_name(patches[i]))
            assert os.path.realpath(link) == os.path.realpath(os.path.dirname(p)), f"link {link} wrong"
        assert seen == {(l, u) for l in (4, 6, 8) for u in (1.0, 4.0)}, f"combinations wrong: {sorted(seen)}"
    with tempfile.TemporaryDirectory() as d:
        # ham_name given as a scalar: the sets in VAR_Hubbard must still be expanded.
        tmpl = os.path.join(d, "parameters")
        open(tmpl, "w").write(SELF_TEST_TEMPLATE.replace('ham_name = {"Hubbard"}', 'ham_name = "Hubbard"'))
        assert len(main(tmpl, RULES, outdir=d)) == 6, "scalar ham_name: VAR_Hubbard sets not expanded"
    print("prepare_sim.py self-test: OK (6 parameter sets, linked lists, rules, symlinks, namelists parse)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-t", "--template", default="parameters", help="template parameters file (default: ./parameters)")
    ap.add_argument("-o", "--outdir", default=".", help="where Temp_<i>/ and the links are created (default: .)")
    ap.add_argument("--no-links", action="store_true", help="do not create the descriptive symbolic links")
    ap.add_argument("--self-test", action="store_true", help="run the generator on a built-in example in a temporary "
                    "directory and check the result")
    args = ap.parse_args()
    if args.self_test:
        self_test()
    else:
        patches = main(args.template, RULES, outdir=args.outdir, links=not args.no_links)
        nper = f90nml.read(args.template).get("var_temp", {}).get("mpi_per_parameter_set")
        print(f"{len(patches)} parameter sets written to {os.path.join(args.outdir, 'Temp_<i>')}")
        if isinstance(nper, int):
            print(f"mpi_per_parameter_set = {nper}: run with mpiexec -n {len(patches) * nper}")
