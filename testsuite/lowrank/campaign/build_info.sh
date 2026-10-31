#!/bin/bash
# build_info.sh FROZEN_BINARY [SOURCE_BINARY]
# Build provenance, run right after the build IN THE SAME SHELL where configure.sh was sourced (it reads ALF_FC,
# ALF_FLAGS_PROG, ALF_LIB, ALF_DIR from the environment), e.g.
#   source configure.sh SUPERMUC-NG MPI NO-INTERACTIVE && make program && testsuite/lowrank/campaign/build_info.sh $WORK/bin/ALF_<tag>.out
# Copies SOURCE_BINARY (default $ALF_DIR/Prog/ALF.out) to FROZEN_BINARY (refuses to overwrite a different existing file)
# and writes FROZEN_BINARY.build_info.txt: KEY=VALUE lines (copied into every job's provenance with prefix build_),
# then raw sections (### ...). Every job verifies the binary's sha256 against binary_sha256 here. bash 3.2, no python.
set -o pipefail
[ $# -ge 1 ] || { sed -n 2,9p "$0" >&2; exit 2; }
dst="$1"
src="${2:-${ALF_DIR:-}/Prog/ALF.out}"
for v in ALF_FC ALF_FLAGS_PROG ALF_LIB ALF_DIR; do
    eval "x=\${$v:-}"
    [ -n "$x" ] || { echo "build_info.sh: $v is not set: source configure.sh in this shell first" >&2; exit 2; }
done
[ -x "$src" ] || { echo "build_info.sh: source binary $src not found" >&2; exit 2; }
. "$(cd "$(dirname "$0")" && pwd)/provenance.sh"

mkdir -p "$(dirname "$dst")"
if [ -e "$dst" ]; then
    if ! cmp -s "$src" "$dst"; then
        echo "build_info.sh: $dst exists and differs from $src; frozen binaries are never overwritten (choose a new name)" >&2
        exit 1
    fi
else
    cp "$src" "$dst" && chmod a-w "$dst" || exit 1
fi
out="$dst.build_info.txt"
[ ! -e "$out" ] || chmod u+w "$out"
: > "$out"
fcver="$($ALF_FC --version 2>&1 | grep -v '^[[:space:]]*$' | head -1)"
if type module >/dev/null 2>&1; then mods="$(module list 2>&1)"; else mods="unavailable(no module command)"; fi
top="$(git -C "$ALF_DIR" rev-parse --show-toplevel 2>/dev/null)"
if [ -n "$top" ] && [ "$(cd "$top" && pwd -P)" = "$(cd "$ALF_DIR" && pwd -P)" ]; then
    commit="$(git -C "$ALF_DIR" rev-parse HEAD)"
    dirty="$(git -C "$ALF_DIR" status --porcelain --untracked-files=no | wc -l | tr -d ' ') modified tracked files"
else
    commit="unavailable(ALF_DIR is not the top of a git tree, e.g. a git archive copy; info has the commit)"; dirty="$commit"
fi
if command -v ldd >/dev/null 2>&1; then libs="$(ldd "$dst" 2>&1)"; libtool=ldd
elif command -v otool >/dev/null 2>&1; then libs="$(otool -L "$dst" 2>&1)"; libtool="otool -L"
else libs="unavailable(no ldd/otool)"; libtool=none; fi
{
    echo "date=$(date '+%Y-%m-%dT%H:%M:%S%z')"
    echo "host=$(hostname)"
    echo "alf_dir=$ALF_DIR"
    echo "source_binary=$src"
    echo "binary=$dst"
    echo "binary_sha256=$(prov_sha256 "$dst")"
    echo "git_commit=$commit"
    echo "git_dirty=$dirty"
    echo "compiler=$ALF_FC"
    echo "compiler_version=$(prov_flat "$fcver")"
    echo "flags_prog=$(prov_flat "$ALF_FLAGS_PROG")"
    echo "alf_lib=$(prov_flat "$ALF_LIB")"
    echo "flags_ext=${ALF_FLAGS_EXT:-none}"
    echo "modules=$(prov_flat "$mods")"
    echo "libs_tool=$libtool"
    echo "libs_mkl=$(printf '%s\n' "$libs" | grep -i -c mkl) lines mention mkl"
    echo "### $ALF_FC --version"
    $ALF_FC --version 2>&1
    echo "### module list"
    printf '%s\n' "$mods"
    echo "### $libtool"
    printf '%s\n' "$libs"
} >> "$out"
chmod a-w "$out"
echo "build_info.sh: $dst (sha256 $(prov_sha256 "$dst")), provenance in $out"
