# TROLL website (pkgdown)

The website <https://troll-model.github.io/TROLL/> is built with
[pkgdown](https://pkgdown.r-lib.org/), the tool used for the rcontroll website, from the
documents of this repository. TROLL is not an R package: `DESCRIPTION` only provides the
metadata pkgdown needs.

| File | Role |
| --- | --- |
| `DESCRIPTION` | site metadata (title, version, authors) |
| `_pkgdown.yml` | navigation bar, article list, home page links |
| `prepare.py` | turns README.md, NEWS.md, docs/*.md, CONTRIBUTING.md, TERMS_OF_REFERENCE.md and sanity/README.md into pkgdown pages, rewriting links between them |
| `build.sh` | runs `prepare.py` and pkgdown; `--api build/docs/html` adds the Doxygen API reference |

Edit the source documents, never the generated files (`index.md`, `NEWS.md`, `vignettes/`,
`docs/`, all ignored by git).

Build locally (requires R with the pkgdown and rmarkdown packages, and Doxygen for the API):

```bash
cmake -S . -B build && cmake --build build --target docs
site/build.sh --api build/docs/html
open site/docs/index.html
```

The workflow `.github/workflows/site.yml` builds the site for every pull request and publishes
it to GitHub Pages on each push to `dev` and each version tag. It needs a one-time setting in
the repository: Settings > Pages > Build and deployment > Source: **GitHub Actions**.
