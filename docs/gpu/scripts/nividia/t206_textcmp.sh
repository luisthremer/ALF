#!/bin/bash
# Ticket 206, nividia: GPU=none identity check of two binaries (HEAD vs working tree), GNU MPI HDF5 build.
#   usage: t206_textcmp.sh [A B]   (default: bin/ALF_gnu_head.out bin/ALF_gnu_new.out, built with t206_build.sh src_head/src_new)
# Compares section sizes, the disassembly of .text (addresses kept: identical addresses mean identical code size) and
# prints the differing lines; differences that are only the constant of a __LINE__ (a "mov $imm,%reg" feeding
# Terminate_on_error) are expected where a source file gained lines (guarded GPU code is preprocessed away, but shifts __LINE__).
W=$HOME/alf-gpu/work/t206
A=${1:-$W/bin/ALF_gnu_head.out}; B=${2:-$W/bin/ALF_gnu_new.out}
. ~/alf-gpu/env/gnu.sh
echo "--- size"; size "$A" "$B"
objdump -d --no-show-raw-insn "$A" | sed 1,3d > "$W/tmp/text_a.txt"
objdump -d --no-show-raw-insn "$B" | sed 1,3d > "$W/tmp/text_b.txt"
echo "--- text lines: $(wc -l < "$W/tmp/text_a.txt") vs $(wc -l < "$W/tmp/text_b.txt")"
diff "$W/tmp/text_a.txt" "$W/tmp/text_b.txt" > "$W/tmp/text_diff.txt"
echo "--- differing lines: $(grep -c '^[<>]' "$W/tmp/text_diff.txt")"
echo "--- differing lines not of the form 'mov/lea $imm' (should be 0):"
grep '^[<>]' "$W/tmp/text_diff.txt" | grep -v -E 'mov +\$0x[0-9a-f]+,|mov +\$0x[0-9a-f]+,%[a-z0-9]+|movl +\$0x[0-9a-f]+,|lea +0x' | head -20
echo "--- sample of differing lines"; grep '^[<>]' "$W/tmp/text_diff.txt" | head -12
echo "--- .rodata/.data differ?"
for s in .rodata .data; do
  objcopy -O binary -j $s "$A" "$W/tmp/sec_a.bin" 2>/dev/null; objcopy -O binary -j $s "$B" "$W/tmp/sec_b.bin" 2>/dev/null
  printf '%s: ' $s; cmp "$W/tmp/sec_a.bin" "$W/tmp/sec_b.bin" > /dev/null 2>&1 && echo identical || echo "differ ($(cmp -l "$W/tmp/sec_a.bin" "$W/tmp/sec_b.bin" 2>/dev/null | wc -l) bytes)"
done
