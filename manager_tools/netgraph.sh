#!/bin/sh
# Open the network graph of this repo and keep it current while this script runs.
# Upstream branches (origin/*) are left out. Ctrl-C stops the updates; the page stays.
here="$(cd "$(dirname "$0")" && pwd)"
python3 "$here/netgraph.py" build --repo "$here/.." --hide-remote origin "$@" || exit 1
open "$here/out/index.html"
exec python3 "$here/netgraph.py" watch --repo "$here/.." --hide-remote origin "$@"
