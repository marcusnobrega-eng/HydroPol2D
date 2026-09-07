#!/bin/zsh

output_dir="/Users/mngomes/Documents/GitHub/HydroPol2D/Applications/Pune_Voronoi_Storm/Outputs/100mm_1h_12h_river_only_d8_verified"
terminal_log="$output_dir/pune-100mm_1h_12h_river_only_d8_verified-terminal.log"
mkdir -p "$output_dir"

{
  echo "START $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  /usr/bin/caffeinate -dimsu /Applications/MATLAB_R2025b.app/bin/matlab -batch \
    "cd('/Users/mngomes/Documents/GitHub/HydroPol2D/Applications/Pune_Voronoi_Storm'); run_pune_voronoi_100mm_12h('pune_river_urban_rural_mesh_d8_2km2_verified','100mm_1h_12h_river_only_d8_verified',[],0.05);"
  exit_code=$?
  echo "FINISH $(date -u '+%Y-%m-%dT%H:%M:%SZ') EXIT_CODE=$exit_code"
} >> "$terminal_log" 2>&1

exit $exit_code
