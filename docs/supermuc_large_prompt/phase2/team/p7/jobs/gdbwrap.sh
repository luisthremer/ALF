#!/bin/bash
export INTELGT_AUTO_ATTACH_DISABLE=1
exec gdb-oneapi -batch -ex "set pagination off" -ex run -ex "bt 30" -ex "info sharedlibrary" --args "$@"
