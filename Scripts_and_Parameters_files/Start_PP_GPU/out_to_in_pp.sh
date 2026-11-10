#!/bin/bash
# Restart helper for PARALLEL_PARAMS runs: in every Temp_<i>/ rename confout_<k>[.h5] to confin_<k>[.h5] (like Start/out_to_in.sh).
# The run then continues exactly (field configuration and random number state are in the confout files) and appends bins to Temp_<i>/data.h5.
shopt -s extglob nullglob
for d in Temp_+([0-9]); do
  ( cd "$d" || exit 1
    if [ -e confout_0.h5 ] && [ -e confout_0 ]; then >&2 echo "$d: both plain text and HDF5 confout present"; exit 1; fi
    for f in confout_+([0-9]).h5 confout_+([0-9]); do mv "$f" "${f/out/in}"; done )
done
