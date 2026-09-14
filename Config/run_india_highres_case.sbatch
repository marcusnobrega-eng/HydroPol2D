#!/bin/bash
# Generic launcher. Submit through submit_india_highres_case.py so memory is
# selected from the case's measured/estimated active-cell count.
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
# Sherlock adds the MATLAB and scratch licenses when it detects the module/use.
# Request only the additional Oak filesystem used by these cases.
#SBATCH --licenses=fs_oak

set -euo pipefail
: "${HP2D_CASE_DIR:?Set HP2D_CASE_DIR}"
: "${HP2D_MODEL_ROOT:?Set HP2D_MODEL_ROOT}"

module --force purge
module load devel math matlab/R2022b python/3.12.1
module load system ffmpeg/7.1

PYTHON=${HP2D_PYTHON:-python3}
readarray -t CASE_VALUES < <("$PYTHON" - "$HP2D_CASE_DIR/case_manifest.json" <<'PY'
import json, sys
from datetime import datetime, timedelta
manifest=json.load(open(sys.argv[1]))
start=datetime.fromisoformat(manifest['simulation_start'])
end=datetime.fromisoformat(manifest['simulation_end'])
smoke=float(__import__('os').environ.get('HP2D_SMOKE_DAYS','0'))
if smoke > 0:
    end=min(end, start+timedelta(days=smoke))
print(start.isoformat())
print(end.isoformat())
print(start.year)
print(manifest['resolution_m'])
PY
)

START=${CASE_VALUES[0]}
END=${CASE_VALUES[1]}
YEAR=${CASE_VALUES[2]}
RESOLUTION=${CASE_VALUES[3]}

"$PYTHON" "$HP2D_MODEL_ROOT/Config/preflight_india_highres_case.py" \
  --case "$HP2D_CASE_DIR" --model-root "$HP2D_MODEL_ROOT" ${HP2D_ALLOW_DIRTY_MODEL:+--allow-dirty-model}

export HYDROPOL_INDIA_CASE_ROOT="$HP2D_CASE_DIR"
export HYDROPOL_INDIA_YEAR="$YEAR"
export HYDROPOL_INDIA_START="$START"
export HYDROPOL_INDIA_END="$END"
export HYDROPOL_INDIA_RESOLUTION_M="$RESOLUTION"
export HYDROPOL_DISABLE_ERA5_DAILY_NETCDF=1
export HYDROPOL_ENABLE_GLOFAS_INFLOW=0
export HYDROPOL_ENABLE_WARMUP=0
export HYDROPOL_USE_GPU=1
export HYDROPOL_DISABLE_GW_2D=0
export HYDROPOL_GW_PERIMETER_MODE=no_flow
export HYDROPOL_RECORD_MAPS_MIN=${HP2D_RECORD_MAPS_MIN:-1440}
if "$PYTHON" -c 'import sys; raise SystemExit(0 if float(sys.argv[1]) > 0 else 1)' "$HP2D_SMOKE_DAYS"; then
    export HYDROPOL_MAX_TIMESTEP_SECONDS=${HP2D_MAX_TIMESTEP_SECONDS:-60}
    if [[ "${HP2D_POSTPROCESS_SMOKE:-0}" == "1" ]]; then
        unset HYDROPOL_SKIP_POSTPROCESS
    else
        export HYDROPOL_SKIP_POSTPROCESS=1
    fi
else
    export HYDROPOL_MAX_TIMESTEP_SECONDS=${HP2D_MAX_TIMESTEP_SECONDS:-300}
    unset HYDROPOL_SKIP_POSTPROCESS
fi
unset HYDROPOL_SYNTHETIC_STRESS_TEST HYDROPOL_SMOKE_DAYS HYDROPOL_DISABLE_GW_LATERAL

echo "Run started: $(date --iso-8601=seconds)"
echo "Case root: $HP2D_CASE_DIR"
echo "Model root: $HP2D_MODEL_ROOT"
echo "Period: $START to $END"
echo "Resolution: $RESOLUTION m"
echo "Memory request: ${SLURM_MEM_PER_NODE:-unknown} MB"

cd "$HP2D_MODEL_ROOT"
export HYDROPOL_PREFLIGHT_ONLY=1
matlab -batch "run(fullfile(pwd,'HydroPol2D_V115.m'))"
export HYDROPOL_PREFLIGHT_ONLY=0
matlab -batch "run(fullfile(pwd,'HydroPol2D_V115.m'))"

echo "Run completed: $(date --iso-8601=seconds)"
