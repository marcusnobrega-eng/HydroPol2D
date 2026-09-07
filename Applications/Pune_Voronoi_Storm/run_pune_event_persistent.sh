#!/bin/zsh
set -eu

case_dir="${0:A:h}"
repository_root="${case_dir:h:h}"
cd "$repository_root"

caffeinate -dimsu /Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "addpath('Applications/Pune_Voronoi_Storm'); run_pune_voronoi_100mm_12h;"
python3 "$case_dir/postprocess_pune_voronoi_event.py"
