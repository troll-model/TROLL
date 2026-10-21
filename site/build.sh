#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Build the TROLL website with pkgdown, from the documents of the repository.
#
# Usage: site/build.sh [--api DOXYGEN_HTML_DIR]
#   --api DIR   Doxygen HTML folder (e.g. build/docs/html): copied into the site as api/, and the
#               XML next to it (build/docs/xml) turned into API reference pages
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
rm -f "$SITE"/vignettes/api-*.Rmd
if [[ -n "$API" && -f "$(dirname "$API")/xml/index.xml" ]]; then
    python3 "$SITE/doxygen_api.py" "$(dirname "$API")/xml" "$SITE"      # API pages from the Doxygen XML
else
    python3 - "$SITE" <<'PY'
import os, sys
site = sys.argv[1]
os.makedirs(os.path.join(site, "vignettes"), exist_ok=True)
for name, title in [("api-index", "API reference"), ("api-tree", "API: Tree"), ("api-species", "API: Species"),
                    ("api-context", "API: Context (model state)"), ("api-functions", "API: functions"),
                    ("api-constants", "API: constants, options and globals"), ("api-types", "API: small types")]:
    open(os.path.join(site, "vignettes", name + ".Rmd"), "w").write(
        f'---\ntitle: "{title}"\n---\n\nThe API pages are generated from the Doxygen XML: build the site with '
        '`site/build.sh --api build/docs/html` after `cmake --build build --target docs`.\n')
PY
fi
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
