# Rollback: restoring the state before the two-lane reorganization (2026-10-29)

Before the verification standard and the two lanes were introduced, everything below was saved. Each section gives
the exact restore commands. Paths under `~/alf-gpu/` are on the test machines; laptop paths are relative to the
directory that holds the ALF clone.

## 1. Git

| Item | Value |
|---|---|
| Tag | `pre-lanes-2026-10-29` (annotated) → `64ecde72` on `gpu-port`, pushed to `fork` |
| 207 WIP | branch `gpu-port-207-wip` on `fork` (on hold, not accepted) |
| Process change | branch `gpu-port-process-lanes` on `fork` (not merged until the owner approves) |

Branch names use `-` instead of `/`: git cannot hold `gpu-port/<x>` next to an existing branch `gpu-port`.

Restore `gpu-port` to the tagged state (only if the process branch was merged and must be undone):

```sh
git fetch fork --tags
git switch gpu-port
git reset --hard pre-lanes-2026-10-29      # local branch
git push --force-with-lease fork gpu-port  # only with the owner's OK: rewrites the public branch
```

Without rewriting history (preferred): `git revert -m 1 <merge-commit>` on `gpu-port`, then push.

Recover the 207 WIP state: `git switch -c gpu-port-207-wip fork/gpu-port-207-wip`.

## 2. Agent definitions and Claude settings (laptop, outside the repo)

`ALF_backups/2026-10-29/` holds:
- `claude_agents/`: copy of `.claude/agents/` at the tag (also recoverable from git: `git checkout pre-lanes-2026-10-29 -- .claude/agents`);
- `claude_memory.tar.gz`: the session memory directory (never committed; stays on the laptop);
- `user_settings.tar.gz`: the user-level Claude settings file; the repo has no `.claude/settings*` file;
- `MANIFEST.sha256`.

Restore the agent definitions:

```sh
git checkout pre-lanes-2026-10-29 -- .claude/agents
```

Restore memory (the directory name is the one in `ALF_backups/2026-10-29/`):

```sh
tar -xzf ALF_backups/2026-10-29/claude_memory.tar.gz -C <claude-project-dir>
```

Agents spawned after a restore use the restored definitions; resumed agents keep their old instructions.

## 3. Test machines

On each machine `~/alf-gpu/backup/env_2026-10-29.tar.gz` holds `~/alf-gpu/env/` and the helper scripts in
`~/alf-gpu/` (and `~/alf-gpu/work/*.sh` on nividia). On `intel` it also holds `backup/oneapi_state_2026-10-29.txt`
(installed oneAPI compiler/MKL directories and the `dpkg` state).

```sh
cd ~/alf-gpu && sha256sum backup/env_2026-10-29.tar.gz   # compare with the value in the gate report
tar -tzf backup/env_2026-10-29.tar.gz                    # inspect first
tar -xzf backup/env_2026-10-29.tar.gz                    # restores env/ and the scripts in place
```

The oneAPI installation itself is not ours to restore (system packages, root).

## 4. 207 evidence archive

`nividia:~/alf-gpu/archive/t207_2026-10-29/` holds a copy of everything left of ticket 207 on 2026-10-29 (logs,
run directories, binaries, source trees from `git archive`, scripts, the raw copies of the paired runs) with
`MANIFEST.sha256`, `bin_sha256.txt` and `tree_binaries_sha256.txt`.

```sh
cd ~/alf-gpu/archive/t207_2026-10-29 && sha256sum -c --quiet MANIFEST.sha256   # no output = intact
```

The originals stay in `~/alf-gpu/work/t207/` until the owner approves their removal.

## 5. NVIDIA MPS (if a daemon was started)

The MPS daemon runs as the ordinary user with its pipe and log directories under `~/alf-gpu/work/mps_trial/`.

```sh
export CUDA_MPS_PIPE_DIRECTORY=~/alf-gpu/work/mps_trial/pipe
export CUDA_MPS_LOG_DIRECTORY=~/alf-gpu/work/mps_trial/log
echo quit | nvidia-cuda-mps-control
pgrep -u "$USER" -f nvidia-cuda-mps || echo "no MPS process left"
nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader   # expect nothing of ours
```

If a server is stuck in a FAULT state: run the `quit` command again with the same pipe directory; if a server process
remains, stop only that process by its PID (it belongs to us), then check `nvidia-smi` again. Never touch other users'
processes.

## 6. Check after any restore

```sh
git describe --tags --exact-match gpu-port   # pre-lanes-2026-10-29 after a hard reset
git diff pre-lanes-2026-10-29 -- CLAUDE.md .claude/agents docs/gpu/tickets/000-template.md   # empty after a restore
```
