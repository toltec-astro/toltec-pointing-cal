#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/calibrate_flux.py"
BACKEND="alma"
DB_PATH=""
TABLE_NAME="flux_estimates"
ARRAY_MAP="0:280,1:220,2:150"
SMA_ALPHA="-0.7"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Loop over observation-number directories and run flux calibration in each one.
Expected ECSV path per obsnum:
  {obsnum}/raw/ppt_commissioning_pointing_{obsnum}_citlali.ecsv

Options:
  -r, --root DIR                 Root directory containing obsnum folders (default: current directory)
  -p, --python-script PATH       Path to calibrate_flux.py (default: ./calibrate_flux.py)
  -b, --backend alma|sqlite      Backend passed to calibrate_flux.py (default: alma)
  -d, --db PATH                  Optional SQLite DB path (required if backend=sqlite)
  -t, --table NAME               SQLite table name (default: flux_estimates)
  -m, --array-map MAP            Array map, e.g. '0:280,1:220,2:150'
      --sma-spectral-index A     Spectral index alpha for SMA fallback (default: -0.7)
  -n, --dry-run                  Print commands without executing
  -h, --help                     Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--root)
      ROOT_DIR="$2"
      shift 2
      ;;
    -p|--python-script)
      PYTHON_SCRIPT="$2"
      shift 2
      ;;
    -b|--backend)
      BACKEND="$2"
      shift 2
      ;;
    -d|--db)
      DB_PATH="$2"
      shift 2
      ;;
    -t|--table)
      TABLE_NAME="$2"
      shift 2
      ;;
    -m|--array-map)
      ARRAY_MAP="$2"
      shift 2
      ;;
    --sma-spectral-index)
      SMA_ALPHA="$2"
      shift 2
      ;;
    -n|--dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$BACKEND" != "alma" && "$BACKEND" != "sqlite" ]]; then
  echo "Error: --backend must be 'alma' or 'sqlite'." >&2
  exit 2
fi

if [[ "$BACKEND" == "sqlite" && -z "$DB_PATH" ]]; then
  echo "Error: --db is required when --backend sqlite is used." >&2
  exit 2
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Error: Python script not found: $PYTHON_SCRIPT" >&2
  exit 2
fi

ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"
if [[ -n "$DB_PATH" ]]; then
  DB_PATH="$(cd "$(dirname "$DB_PATH")" && pwd)/$(basename "$DB_PATH")"
fi

shopt -s nullglob
processed=0

for obs_dir in "$ROOT_DIR"/[0-9]*; do
  [[ -d "$obs_dir" ]] || continue

  obsnum="$(basename "$obs_dir")"
  [[ "$obsnum" =~ ^[0-9]+$ ]] || continue

  ecsv_rel="raw/ppt_commissioning_pointing_${obsnum}_citlali.ecsv"
  ecsv_path="$obs_dir/$ecsv_rel"

  if [[ ! -f "$ecsv_path" ]]; then
    continue
  fi

  output_file="cal_ratios_${obsnum}.json"

  cmd=(python "$PYTHON_SCRIPT" "$ecsv_rel")
  if [[ -n "$DB_PATH" ]]; then
    cmd+=("$DB_PATH")
  fi

  cmd+=(
    --backend "$BACKEND"
    --table "$TABLE_NAME"
    --array-map "$ARRAY_MAP"
    --sma-spectral-index "$SMA_ALPHA"
    --output "$output_file"
  )

  echo "[obsnum=$obsnum] running calibration in $obs_dir"
  if [[ $DRY_RUN -eq 1 ]]; then
    (cd "$obs_dir" && printf '  DRY RUN: '; printf '%q ' "${cmd[@]}"; printf '\n')
  else
    (cd "$obs_dir" && "${cmd[@]}")
  fi

  processed=$((processed + 1))
done

if [[ $processed -eq 0 ]]; then
  echo "No matching obsnum directories with expected ECSV files were found under: $ROOT_DIR" >&2
  exit 1
fi

echo "Finished. Processed $processed observation director$( [[ $processed -eq 1 ]] && echo 'y' || echo 'ies' )."
