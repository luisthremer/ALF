# Identical helper scripts (checked 2026-11-13)

These files are byte-identical copies (same `shasum`). Nothing was removed; the first path of
each line is the copy to keep, the others are removal candidates (owner decision).
Re-check with:
`find docs/gpu -type f \( -name '*.sh' -o -name '*.py' \) -exec shasum {} + | sort`

| Keep | Identical copies |
|---|---|
| `docs/gpu/scripts/nividia/t207_firstdiv.py` | `scripts/t207_final/`, `scripts/t208_1/`, `scripts/t208_2/`, `scripts/t208_3/` (same name) |
| `docs/gpu/scripts/defect/compare.py` | `sessions/supermuc_large_prompt/jobs/compare.py`, `.../phase2/team/compare.py` |
| `docs/gpu/sessions/supermuc_large_prompt/jobs/env_intel_pvc.sh` | `.../diag/env_intel_pvc.sh`, `.../phase2/team/env_intel_pvc.sh` |
| `docs/gpu/sessions/supermuc_large_prompt/jobs/bind_gpu_per_card.sh` | `.../phase2/team/bind_gpu_per_card.sh` |
| `docs/gpu/sessions/supermuc_large_prompt/jobs/run_ref_sng.py` | `.../phase2/team/run_ref_sng.py` |
| `docs/gpu/sessions/supermuc_large_prompt/phase2/team/comb/build_comb.sh` | `.../phase2/team/p5/build_comb.sh` |
| `docs/gpu/scripts/nividia/t207_d5_h5cmp.sh` | `scripts/nividia/t207_1/t207_1_d5_h5cmp.sh` |
| `docs/gpu/scripts/nividia/t207_h5cmp.sh` | `scripts/nividia/t205_h5cmp.sh` |

Scripts called from job files or other scripts by relative path may depend on their copy;
check with `grep -rn <name>` before removing one.
