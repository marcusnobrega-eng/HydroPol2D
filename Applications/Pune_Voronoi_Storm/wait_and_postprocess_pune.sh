#!/bin/zsh
set -eu

matlab_pid="$1"
case_dir="${0:A:h}"
result_file="$case_dir/Outputs/100mm_1h_12h/pune-100mm-1h-12h-results.nc"

while kill -0 "$matlab_pid" 2>/dev/null; do
  sleep 30
done

if [[ ! -f "$result_file" ]]; then
  print -u2 "MATLAB ended without creating $result_file"
  exit 1
fi

python3 "$case_dir/postprocess_pune_voronoi_event.py"
