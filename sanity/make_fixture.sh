#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Regenerate sanity/fixtures/stand_1ha: a 1 ha (100 x 100 m) stand grown from
# bare ground with the ./example inputs, saved as an inventory (-f) plus soil
# water content (-w) that TROLL can restart from.
#
# Usage: sanity/make_fixture.sh --troll PATH [--years N] [--chunk-years N]
#
# The run is split into chunks that restart from the previous chunk's
# final_pattern / final_SWC3D. This keeps each process short and exercises the
# inventory-restart path.
#
# The committed fixture was made with --years 15 (712 trees > 10 cm/ha,
# AGB ~88 t/ha). A longer spin-up (e.g. --years 300) gives a closer-to-
# equilibrium stand that also exercises large-tree dynamics; regenerate it on
# a workstation and commit it together with an updated reference if desired.
# ----------------------------------------------------------------------------
set -euo pipefail

SANITY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SANITY_DIR/.." && pwd)"
EX="$ROOT_DIR/example"
OUT="$SANITY_DIR/fixtures/stand_1ha"

TROLL_EXE=""
YEARS=15
CHUNK=5
while [[ $# -gt 0 ]]; do
    case "$1" in
    --troll) TROLL_EXE="$2"; shift 2 ;;
    --years) YEARS="$2"; shift 2 ;;
    --chunk-years) CHUNK="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -x "$TROLL_EXE" ]] || { echo "--troll PATH (an executable) is required" >&2; exit 2; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

sed -E -e 's/^cols[[:space:]]+[0-9]+/cols\t100/' -e 's/^rows[[:space:]]+[0-9]+/rows\t100/' \
    -e "s/^nbiter[[:space:]]+[0-9]+/nbiter\t$((CHUNK * 365))/" -e 's/^_NONRANDOM[[:space:]]+[0-9]+/_NONRANDOM\t1/' \
    "$EX/global_inputs.txt" > "$WORK/global_inputs.txt"

done_years=0
restart=()
chunk=0
while [[ $done_years -lt $YEARS ]]; do
    chunk=$((chunk + 1))
    echo "Chunk $chunk: years $done_years -> $((done_years + CHUNK))"
    "$TROLL_EXE" -i"$WORK/global_inputs.txt" -s"$EX/species.txt" -m"$EX/daily_climate.txt" \
        -d"$EX/halfhourly_climate.txt" -p"$EX/soil.txt" -o"$WORK/c$chunk" ${restart[@]+"${restart[@]}"} > "$WORK/c$chunk.log" 2>&1
    grep "Trees (1/ha)" "$WORK/c$chunk.log" | tail -1 | cut -c1-120
    restart=("-f$WORK/c${chunk}_0_final_pattern.txt" "-w$WORK/c${chunk}_0_final_SWC3D.txt")
    done_years=$((done_years + CHUNK))
done

mkdir -p "$OUT"
gzip -9 -c "$WORK/c${chunk}_0_final_pattern.txt" > "$OUT/forest.txt.gz"
cp "$WORK/c${chunk}_0_final_SWC3D.txt" "$OUT/SWC3D.txt"
sed -E 's/^nbiter[[:space:]]+[0-9]+/nbiter\t365/' "$WORK/global_inputs.txt" > "$OUT/global_inputs.txt"
echo "Fixture written to $OUT (stand age: $done_years years)."
echo "Reference outputs for --case stand are now stale: rerun the sanity check to rebuild the cache."
