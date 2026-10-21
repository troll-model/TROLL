#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Build the TROLL website with pkgdown, from the documents of the repository.
#
# Usage: site/build.sh [--api DOXYGEN_HTML_DIR]
#   --api DIR   copy the Doxygen HTML (e.g. build/docs/html) into the site as api/
#
# The pages are generated from README.md, docs/*.md, NEWS.md, CONTRIBUTING.md,
# TERMS_OF_REFERENCE.md and sanity/README.md, so each text exists only once.
# Output: site/docs/ (open site/docs/index.html). Requires R with pkgdown.
# pkgdown.internet = FALSE skips the CRAN lookup (TROLL is not an R package on CRAN).
# ----------------------------------------------------------------------------
set -euo pipefail
SITE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SITE/.." && pwd)"
API=""
while [[ $# -gt 0 ]]; do
    case "$1" in
    --api) API="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

python3 "$SITE/prepare.py" "$ROOT" "$SITE"
Rscript -e 'options(pkgdown.internet = FALSE); pkgdown::build_site(pkg = commandArgs(TRUE)[1], preview = FALSE, install = FALSE, examples = FALSE, new_process = FALSE)' "$SITE"

if [[ -n "$API" ]]; then
    rm -rf "$SITE/docs/api"
    cp -R "$API" "$SITE/docs/api"
    echo "API reference copied to site/docs/api"
else
    mkdir -p "$SITE/docs/api"
    cat > "$SITE/docs/api/index.html" <<'HTML'
<!doctype html><meta charset="utf-8"><title>TROLL API reference</title>
<p>The C++ API reference is generated with Doxygen: run <code>cmake --build build --target docs</code>
and open <code>build/docs/html/index.html</code>, or build the site with <code>site/build.sh --api build/docs/html</code>.</p>
HTML
fi
echo "Website built in $SITE/docs (open $SITE/docs/index.html)"
