import f90nml
import itertools
import os
import subprocess
from collections import defaultdict
import json

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

def unflatten(data, raw_combos, conditions):
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

def main(original_param_path, rules):
    original_nml = f90nml.read(original_param_path)
    extracted_information = cleanup_dirty_list(original_nml)
    data = extracted_information

    flattend = flatten_dict(data)
    grouped = group_something(flatten_dict(data))
    linked = zip_linked(flattend, grouped)
    l_of_patches = safely_deduplicate(unflatten(data, cartesian_product(linked), rules))
    
    for i, el in enumerate(l_of_patches):
        sim_path = os.path.join(".", f"Temp_{i}")
        param_path = os.path.join(sim_path, "parameters")
        
        subprocess.run(["mkdir", "-p", sim_path], check=True)
        subprocess.run(["cp", "parameters", sim_path], check=True)
        
        apply_patch(original_param_path, el, param_path)

        link_name = f"sim_{i}"
        if os.path.lexists(link_name):
            os.unlink(link_name)
        subprocess.run(["ln", "-s", sim_path, link_name], check=True)


if __name__=="__main__":
    rules = {
        'var_hubbard': ('ham_name', 'Hubbard'),
        'var_tv': ('ham_name', 'tV')
    }

    main("parameters", rules)


    #TODO: def automatic_rules()
    #TODO: link names
    #TODO: MPI PER PARAMETER SET klären!
    #TODO: MODEL to list!
