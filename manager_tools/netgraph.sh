#!/bin/sh
# Open the network graph of this repo and keep it current while this script runs.
# Which branches are shown: netgraph.json (default: the current Code-modularity line, no origin/*). Ctrl-C stops the updates; the page stays.
here="$(cd "$(dirname "$0")" && pwd)"
python3 "$here/netgraph.py" build --repo "$here/.." "$@" || exit 1
open "$here/out/index.html"
exec python3 "$here/netgraph.py" watch --repo "$here/.." "$@"
