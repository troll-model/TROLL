# TROLL

[![DOI](https://zenodo.org/badge/49142579.svg)](https://doi.org/10.5281/zenodo.14013146)

TROLL is an individual- and trait-based model of tropical forest dynamics. It simulates every
tree in a forest stand of 1 to 100 ha on a 1 m resolution 3D grid, and tracks forest structure,
carbon and water fluxes, and species-level biodiversity in the same simulation.

This repository holds **TROLL 4.1.11**, the current implementation of the **TROLL 4.0 model**
described in Maréchaux et al. (2025) and evaluated in Schmitt et al. (2025). Versions 4.1.x
reorganise, test, document and accelerate the code; two of them change simulation results
slightly (4.1.1: single-precision leaf physiology; 4.1.2: a correction to the light-limited
photosynthesis rate). See [NEWS.md](NEWS.md) for every change.

Website: <https://troll-model.github.io/TROLL/> · R interface: [rcontroll](https://sylvainschmitt.github.io/rcontroll/)

## How to cite

Please cite the model description, and the evaluation if relevant, together with the version
of the code you used ([CITATION.cff](CITATION.cff); GitHub shows it under "Cite this repository"):

- Maréchaux, I., Fischer, F. J., Schmitt, S., & Chave, J. (2025). TROLL 4.0: representing water and carbon fluxes, leaf phenology, and intraspecific trait variation in a mixed-species individual-based forest dynamics model – Part 1: Model description. *Geoscientific Model Development*, 18(16), 5143–5204. https://doi.org/10.5194/gmd-18-5143-2025
- Schmitt, S., Fischer, F. J., Ball, J. G. C., Barbier, N., Boisseaux, M., Bonal, D., Burban, B., Chen, X., Derroire, G., Lichstein, J. W., Nemetschek, D., Restrepo-Coupe, N., Saleska, S., Sellan, G., Verley, P., Vincent, G., Ziegler, C., Chave, J., & Maréchaux, I. (2025). TROLL 4.0: representing water and carbon fluxes, leaf phenology, and intraspecific trait variation in a mixed-species individual-based forest dynamics model – Part 2: Model evaluation for two Amazonian sites. *Geoscientific Model Development*, 18(16), 5205–5243. https://doi.org/10.5194/gmd-18-5205-2025

Earlier versions and the R interface are referenced in the [Terms of Reference](TERMS_OF_REFERENCE.md#5-references).

## Installation

Requirements:

| Tool | Needed for | Ubuntu | macOS (Homebrew) |
| --- | --- | --- | --- |
| C++17 compiler | building | `sudo apt install g++` (or `clang`) | Xcode Command Line Tools: `xcode-select --install` |
| CMake ≥ 3.13 | building | `sudo apt install cmake` | `brew install cmake` |
| GNU Scientific Library | building | `sudo apt install libgsl-dev` | `brew install gsl` |
| OpenMP (optional) | multi-core runs | included with g++ | `brew install libomp` |
| Python 3 | tests | `sudo apt install python3` | included with macOS |
| Doxygen (optional) | API documentation | `sudo apt install doxygen` | `brew install doxygen` |

Build, from the top-level folder of the repository:

```bash
cmake -S . -B build
cmake --build build
```

The executable is `build/TROLL`.

On macOS with Apple Silicon, use the native Homebrew in `/opt/homebrew`; an Intel Homebrew in
`/usr/local` produces x86_64 builds that fail to link with arm64 libraries. Homebrew's `libomp`
must be pointed to explicitly:

```bash
cmake -S . -B build -DOpenMP_ROOT=$(brew --prefix libomp)
```

### Build options

| CMake option | Default | Effect |
| --- | --- | --- |
| `-DTROLL_WARNINGS_AS_ERRORS=ON` | `OFF` | `-Werror` (used in CI) |
| `-DTROLL_EXTRA_WARNINGS=ON` | `OFF` | adds `-Wextra -Wshadow -Wconversion` (advisory) |
| `-DTROLL_SANITIZE="address;undefined"` | empty | builds with sanitizers (use with `-DCMAKE_BUILD_TYPE=Debug`) |
| `-DTROLL_OPENMP=OFF` | `ON` | parallel leaf fluxes and LAImax table (OpenMP); results do not depend on the thread count |
| `-DTROLL_LTO=OFF` | `ON` | link-time optimisation |
| `-DTROLL_DEFINES="..."` | empty | compile-time model options, see below |

Model options are listed in `include/troll_defines.hpp`. Options that are on by default are
switched off with `TROLL_NO_<OPTION>`; options that are off by default are switched on by naming
them, e.g. `cmake -S . -B build -DTROLL_DEFINES="LAI_gradient;CACHE_LUT"`. CI builds each option;
`TROLL_NO_WATER`, `TROLL_NO_CROWN_UMBRELLA`, `TROLL_NO_FULL_CLIMATE`, `Output_ABC` and
`TRACK_INDIVIDUALS` currently do not compile.

## Running TROLL

TROLL reads five input files: global parameters, species parameters, soil, and climate at daily
and half-hourly resolution ([docs/inputs.md](docs/inputs.md)). With the example inputs:

```bash
./build/TROLL -i./example/global_inputs.txt -s./example/species.txt -m./example/daily_climate.txt \
  -d./example/halfhourly_climate.txt -p./example/soil.txt -o./TROLL_test
```

Optional arguments: `-f<file>` and `-w<file>` initialise the forest and soil water from an
inventory (e.g. the `_final_pattern.txt` and `_final_SWC3D.txt` outputs of a previous run),
`-n<index>` sets the run index used in output file names, and `-r<seed>` (or `--seed=<seed>`)
sets the random seed. Without `-r`, the seed is 1 if `_NONRANDOM` is 1 in the global inputs and
time-based otherwise. The TROLL version and the seed are written to the `_info.txt` output so
every run can be reproduced. The outputs are described in [docs/outputs.md](docs/outputs.md).

TROLL uses all available cores by default; set `OMP_NUM_THREADS` to choose the number of threads
(e.g. `OMP_NUM_THREADS=1` when running many simulations side by side). Outputs are
bit-identical whatever the number of threads.

## Tests

```bash
cmake --build build --target check          # unit tests + regression checks (about 2 minutes)
cmake --build build --target threads-check  # outputs must be bit-identical with 1 and 4 threads
cmake --build build --target sanity1000     # 1000 days on an established stand (about 4 minutes)
```

The regression checks compare against reference outputs stored in `sanity/golden/`. On Linux
x86-64 with GCC 13, results reproduce them bit for bit; the default tolerances allow for other
compilers and math libraries. See [sanity/README.md](sanity/README.md).

## Documentation

| Document | Content |
| --- | --- |
| [docs/inputs.md](docs/inputs.md) | input files, columns and units |
| [docs/outputs.md](docs/outputs.md) | output files, columns and units |
| [docs/equations.md](docs/equations.md) | where each equation of the model description is implemented, with constants checked |
| [docs/open-questions.md](docs/open-questions.md) | discrepancies between code, papers and files, for the authors |
| [CONTRIBUTING.md](CONTRIBUTING.md) | build, test and documentation rules for contributors |
| [NEWS.md](NEWS.md) | changes in each version |
| API reference | `cmake --build build --target docs` (requires Doxygen; HTML in `build/docs/html`), or on the website |

## R interface

The R package [rcontroll](https://github.com/sylvainschmitt/rcontroll) prepares inputs, runs
TROLL and analyses its outputs (Schmitt et al. 2023;
[website](https://sylvainschmitt.github.io/rcontroll/)). rcontroll compiles its own copy of the
TROLL code; check which TROLL version it embeds before comparing results with this repository.
The version supporting TROLL 4.0 is installed with:

```R
devtools::install_github("sylvainschmitt/rcontroll", ref = "TROLLV4")
```

## Help and contributions

Please use the GitHub [issues](https://github.com/troll-model/TROLL/issues) for any question
about TROLL. Contributions are welcome: see [CONTRIBUTING.md](CONTRIBUTING.md).

## Core development team

The project is governed by the core development team described in the
[Terms of Reference](TERMS_OF_REFERENCE.md), who can be contacted at dev@troll-model.org:

- **Isabelle Maréchaux**, INRAE, AMAP, Montpellier, France (project lead)
- **Fabian Jörg Fischer**, Technical University of Munich (TUM), Germany
- **Sylvain Schmitt**, CIRAD, Montpellier, France
- **Philippe Verley**, IRD, AMAP, Montpellier, France
- **Jérôme Chave**, CNRS, CRBE, Toulouse, France

## Licence

TROLL is distributed under the GNU General Public License v3.0 ([LICENSE](LICENSE)).
