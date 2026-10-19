#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Time TROLL on the established-stand fixture for several thread counts.
#
# Usage: sanity/bench.sh --troll PATH [--nbiter N] [--threads "1 2 4 8"]
#
# Prints wall time, time in the simulation loop (from _info.txt) and start-up
# + output time per thread count, and checks that outputs are identical to the
# single-thread run.
# ----------------------------------------------------------------------------
set -euo pipefail
SANITY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EX="$SANITY_DIR/../example"
FX="$SANITY_DIR/fixtures/stand_1ha"
TROLL_EXE=""; NBITER=60; THREADS="1 2 4 8"
while [[ $# -gt 0 ]]; do
    case "$1" in
    --troll) TROLL_EXE="$2"; shift 2 ;;
    --nbiter) NBITER="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -x "$TROLL_EXE" ]] || { echo "--troll PATH (an executable) is required" >&2; exit 2; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
gunzip -c "$FX/forest.txt.gz" > "$WORK/forest.txt"
sed -E "s/^nbiter[[:space:]]+[0-9]+/nbiter\t$NBITER/" "$FX/global_inputs.txt" > "$WORK/global_inputs.txt"

printf "%-8s %10s %10s %14s %s\n" threads "wall (s)" "loop (s)" "startup+io (s)" "outputs vs 1 thread"
for t in $THREADS; do
    mkdir -p "$WORK/t$t"
    start=$(date +%s.%N)
    OMP_NUM_THREADS=$t "$TROLL_EXE" -i"$WORK/global_inputs.txt" -s"$EX/species.txt" -m"$EX/daily_climate.txt" \
        -d"$EX/halfhourly_climate.txt" -p"$EX/soil.txt" -f"$WORK/forest.txt" -w"$FX/SWC3D.txt" -o"$WORK/t$t/run" > "$WORK/t$t.log" 2>&1
    wall=$(python3 -c "import time; print(f'{time.time() - $start:.2f}')")
    loop=$(awk '/Average computation time/ {print $5}' "$WORK/t$t/run_0_info.txt")
    rest=$(python3 -c "print(f'{$wall - $loop:.2f}')")
    same="-"
    if [[ -d "$WORK/t1" && "$t" != 1 ]]; then
        python3 "$SANITY_DIR/compare.py" "$WORK/t1" "$WORK/t$t" > "$WORK/cmp$t.txt" && same="identical" || same="DIFFERENT (see compare.py)"
    fi
    printf "%-8s %10s %10s %14s %s\n" "$t" "$wall" "$loop" "$rest" "$same"
done
