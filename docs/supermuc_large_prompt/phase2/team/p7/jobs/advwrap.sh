#!/bin/bash
source /lrz/sys/intel/oneapi_2025.3.1/advisor/latest/advisor-vars.sh >/dev/null 2>&1
exec timeout ${ADV_TIMEOUT:-420} advisor --collect=${ADV_COLLECT:-survey} ${ADV_FLAGS} --profile-gpu --project-dir=${ADV_PROJ:-$PWD/adv} -- "$@"
