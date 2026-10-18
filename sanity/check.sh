#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# TROLL regression ("sanity") check
#
# Runs a pinned reference version of TROLL and the current build on the same
# inputs, then compares every output file (bit-exact by default).
#
# Usage: sanity/check.sh --troll PATH [options]
#   --troll PATH    TROLL executable to test (required; CMake passes it)
#   --case NAME     example : 4 ha from bare ground, ./example inputs (default)
#                   stand   : 1 ha established stand from sanity/fixtures/stand_1ha,
#                             exercising adult growth, competition, mortality, treefall
#   --nbiter N      number of daily timesteps (default 20 for example, 365 for stand)
#   --seed N        RNG seed passed with -rN (needs a reference that supports -r)
#   --ref SHA       reference commit (default: $TROLL_SANITY_REF or the pinned SHA below)
#   --rtol X        relative tolerance for numeric fields (default 0 = bit-exact)
#   --atol X        absolute tolerance for numeric fields (default 0)
#
# The reference build and its outputs are cached in sanity/.ref/<sha>/.
# ----------------------------------------------------------------------------
set -euo pipefail

# Pinned reference: dev at 4de3273 ("41 terms of reference"), i.e. after the
# monolith split, which itself was validated against main (v4.0).
# Update deliberately, in its own commit, when outputs are *meant* to change.
PINNED_REF="4de327370e3dae47f8cf999eec91c4d30487b65b"
REPO_URL="${TROLL_SANITY_REPO:-https://github.com/troll-model/TROLL.git}"

SANITY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SANITY_DIR/.." && pwd)"

TROLL_EXE=""
CASE="example"
NBITER=""
SEED=""
REF="${TROLL_SANITY_REF:-$PINNED_REF}"
RTOL="0"
ATOL="0"

while [[ $# -gt 0 ]]; do
    case "$1" in
    --troll) TROLL_EXE="$2"; shift 2 ;;
    --case) CASE="$2"; shift 2 ;;
    --nbiter) NBITER="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --rtol) RTOL="$2"; shift 2 ;;
    --atol) ATOL="$2"; shift 2 ;;
    -h | --help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; NC='\033[0m'
info() { echo -e "${YELLOW}$*${NC}"; }

[[ -n "$TROLL_EXE" ]] || { echo "--troll PATH is required" >&2; exit 2; }
TROLL_EXE="$(cd "$(dirname "$TROLL_EXE")" && pwd)/$(basename "$TROLL_EXE")"
[[ -x "$TROLL_EXE" ]] || { echo "Not an executable: $TROLL_EXE" >&2; exit 2; }

case "$CASE" in
example) NBITER="${NBITER:-20}" ;;
stand) NBITER="${NBITER:-365}" ;;
*) echo "Unknown case: $CASE (expected example or stand)" >&2; exit 2 ;;
esac

info "=== TROLL SANITY CHECK: case=$CASE nbiter=$NBITER seed=${SEED:-input file} ref=${REF:0:10} ==="

# ----------------------------------------------------------------------------
# 1. Build the reference (cached per commit)
# ----------------------------------------------------------------------------
REF_DIR="$SANITY_DIR/.ref/$REF"
REF_EXE="$REF_DIR/TROLL_ref"
if [[ ! -x "$REF_EXE" ]]; then
    info "Building reference $REF ..."
    rm -rf "$REF_DIR/src" "$REF_DIR/build"
    mkdir -p "$REF_DIR"
    git clone --quiet "$REPO_URL" "$REF_DIR/src"
    git -C "$REF_DIR/src" checkout --quiet "$REF"
    if [[ -f "$REF_DIR/src/CMakeLists.txt" ]]; then
        cmake -S "$REF_DIR/src" -B "$REF_DIR/build" -DCMAKE_BUILD_TYPE=Release > "$REF_DIR/build.log"
        cmake --build "$REF_DIR/build" --target TROLL -j 4 >> "$REF_DIR/build.log"
        cp "$REF_DIR/build/TROLL" "$REF_EXE"
    else # pre-CMake monolith (main / v4.0)
        if [[ "$(uname -s)" == "Darwin" ]]; then GSL_PREFIX="$(brew --prefix gsl)"; else GSL_PREFIX="$(pkg-config --variable=prefix gsl)"; fi
        g++ -w -O3 -o "$REF_EXE" "$REF_DIR/src/mainTROLL4.0.cpp" -I"$GSL_PREFIX/include" -L"$GSL_PREFIX/lib" -lgsl -lgslcblas -lm
    fi
    echo -e "${GREEN}✔ reference built${NC}"
fi

if [[ -n "$SEED" ]] && ! grep -rqs -- "--seed" "$REF_DIR/src/src" "$REF_DIR/src/mainTROLL4.0.cpp"; then
    echo "The reference ${REF:0:10} has no -r/--seed option; use --ref with a newer commit or drop --seed." >&2
    exit 2
fi

# ----------------------------------------------------------------------------
# 2. Prepare inputs for the selected case
# ----------------------------------------------------------------------------
CASE_ID="${CASE}_n${NBITER}_s${SEED:-default}"
WORK="$SANITY_DIR/test_runs/$CASE_ID"
rm -rf "$WORK"
mkdir -p "$WORK/inputs" "$WORK/new"

EX="$ROOT_DIR/example"
EXTRA_ARGS=()
if [[ "$CASE" == "example" ]]; then
    GLOBAL_SRC="$EX/global_inputs.txt"
else
    FX="$SANITY_DIR/fixtures/stand_1ha"
    GLOBAL_SRC="$FX/global_inputs.txt"
    gunzip -c "$FX/forest.txt.gz" > "$WORK/inputs/forest.txt"
    cp "$FX/SWC3D.txt" "$WORK/inputs/SWC3D.txt"
    EXTRA_ARGS+=("-f$WORK/inputs/forest.txt" "-w$WORK/inputs/SWC3D.txt")
fi
sed -E "s/^nbiter[[:space:]]+[0-9]+/nbiter\t$NBITER/" "$GLOBAL_SRC" > "$WORK/inputs/global_inputs.txt"
grep -qE "^nbiter[[:space:]]+$NBITER([[:space:]]|$)" "$WORK/inputs/global_inputs.txt" || { echo "Could not set nbiter" >&2; exit 1; }
if [[ -n "$SEED" ]]; then EXTRA_ARGS+=("-r$SEED"); fi

run_troll() { # $1 executable, $2 output prefix, $3 log file
    "$1" -i"$WORK/inputs/global_inputs.txt" -s"$EX/species.txt" -m"$EX/daily_climate.txt" \
        -d"$EX/halfhourly_climate.txt" -p"$EX/soil.txt" -o"$2" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} > "$3" 2>&1
}

# Reference outputs are cached as well, keyed on everything that defines the run
INPUT_HASH="$(cat "$WORK/inputs/global_inputs.txt" "$EX"/*.txt "$SANITY_DIR/fixtures/stand_1ha"/* | cksum | cut -d' ' -f1)"
REF_RUN="$REF_DIR/runs/${CASE_ID}_$INPUT_HASH"
if [[ ! -f "$REF_RUN/.complete" ]]; then
    info "Running reference ..."
    rm -rf "$REF_RUN"; mkdir -p "$REF_RUN/out"
    run_troll "$REF_EXE" "$REF_RUN/out/test" "$REF_RUN/run.log" || { tail -20 "$REF_RUN/run.log"; exit 1; }
    touch "$REF_RUN/.complete"
fi

# ----------------------------------------------------------------------------
# 3. Run the current build and compare
# ----------------------------------------------------------------------------
info "Running current build ..."
run_troll "$TROLL_EXE" "$WORK/new/test" "$WORK/run.log" || { tail -20 "$WORK/run.log"; exit 1; }
echo -e "${GREEN}✔ run complete${NC}"

info "Comparing outputs ..."
if python3 "$SANITY_DIR/compare.py" "$REF_RUN/out" "$WORK/new" --rtol "$RTOL" --atol "$ATOL"; then
    echo -e "${GREEN}SANITY CHECK PASSED ($CASE_ID)${NC}"
else
    echo -e "${RED}SANITY CHECK FAILED ($CASE_ID)${NC}"
    exit 1
fi
