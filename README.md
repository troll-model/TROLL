# TROLL

[![DOI](https://zenodo.org/badge/49142579.svg)](https://doi.org/10.5281/zenodo.14013146)

TROLL 4.0 is an individual- and trait-based model of forest dynamics. The model and its evaluation are described in two companion papers in Geoscientific Model Development:

- Maréchaux I., Fischer F.J., Schmitt S., Chave J. (2025). TROLL 4.0: representing water and carbon fluxes, leaf phenology, and intraspecific trait variation in a mixed-species individual-based forest dynamics model – Part 1: Model description. *Geosci. Model Dev.*, 18, 5143–5204. https://doi.org/10.5194/gmd-18-5143-2025
- Schmitt S., Fischer F.J., *et al.* (2025). TROLL 4.0 … – Part 2: Model evaluation for two Amazonian sites. *Geosci. Model Dev.*, 18, 5205–5243. https://doi.org/10.5194/gmd-18-5205-2025

## Installation

Required:
- CMake Build System, see https://cmake.org
  - Ubuntu: `sudo apt install cmake`
  - macOS: `brew install cmake`
- GNU Scientific Library (GSL) installation, see https://www.gnu.org/software/gsl/
  - Ubuntu: `sudo apt install libgsl-dev`
  - macOS: `brew install gsl`
- C++ compiler, both [`g++`](https://gcc.gnu.org/) (GNU Compiler Collection) and [`clang++`](https://clang.llvm.org/) (C language family frontend for LLVM) work fine
  - Ubuntu: 
    - Default install g++  `sudo apt install g++`
    - Alternatively install clang++ `sudo apt install clang++`
  - macOS: 
    - Default install Xcode Command Line Tools (includes `clang++`): `xcode-select --install`
    - Alternatively install `g++` via Homebrew `brew install gcc`

The following command line can be used to compile the code:

```bash
cd build/
cmake ..
make
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

Model options are listed in `include/troll_defines.hpp`. Options that are on by default are switched off with `TROLL_NO_<OPTION>`; options that are off by default are switched on by naming them, e.g. `cmake -DTROLL_DEFINES="LAI_gradient;CACHE_LUT" ..`. CI builds each option; `TROLL_NO_WATER`, `TROLL_NO_CROWN_UMBRELLA`, `TROLL_NO_FULL_CLIMATE`, `Output_ABC` and `TRACK_INDIVIDUALS` currently do not compile.

### Tests

```bash
cd build/
make unit-tests    # process-level unit tests (tests/)
make sanity        # 20 days from bare ground, against stored reference outputs
make sanity-stand  # 1 year on an established 1 ha stand
make check         # all three
make sanity1000    # 1000 days on the established stand (~4 min)
make threads-check # outputs must be bit-identical with 1 and 4 threads
make golden        # regenerate the stored reference outputs (only when results are meant to change)
```

The regression checks compare against reference outputs stored in `sanity/golden/`. On Linux x86-64 with GCC 13, results reproduce them bit for bit; the default tolerances (CMake variables `TROLL_SANITY_RTOL_EXAMPLE`, `TROLL_SANITY_RTOL_LONG`, `TROLL_SANITY_ATOL`, `TROLL_SANITY_FILES_LONG`) allow for other compilers and math libraries. See [sanity/README.md](sanity/README.md).

The leaf kernel (`FluxesLeaf`, `Photosyn`, `QUAD`) is single precision and compiled with `-Wdouble-promotion` as an error (`include/kernel_precision.hpp`): use the `*_F` constants of `constants.hpp` and `f`-suffixed literals there. Building with `-ffast-math` is refused.

### Documentation

- [docs/outputs.md](docs/outputs.md): every output file, its columns and units.
- API documentation: `make docs` generates it in `build/docs/html` (requires Doxygen); `make docs-check` also reports coverage and fails on documentation errors.

## Get started

The code of TROLL 4.0 is written in C++ and requires five input files to run a simulation: (i) global parameters (global_inputs), (ii) species parameters (species), (iii) soil characteristics (soil), and finally, meteorological drivers varying at (iv) half-hour (halfhourly_climate) and (v) daily step (daily_climate).

The following command line can be used to run the code:

```bash
./build/TROLL -i./example/global_inputs.txt -s./example/species.txt -m./example/daily_climate.txt -d./example/halfhourly_climate.txt -p./example/soil.txt -o./TROLLv4_test
```

Optional arguments: `-f<file>` and `-w<file>` initialise the forest and soil water from an inventory (e.g. the `_final_pattern.txt` and `_final_SWC3D.txt` outputs of a previous run), `-n<index>` sets the run index used in output file names, and `-r<seed>` (or `--seed=<seed>`) sets the random seed. Without `-r`, the seed is 1 if `_NONRANDOM` is 1 in the global inputs and time-based otherwise. The seed actually used is written to the `_info.txt` output so every run can be reproduced.

TROLL uses all available cores by default; set `OMP_NUM_THREADS` to choose the number of threads (e.g. `OMP_NUM_THREADS=1` when running many simulations side by side). Outputs are bit-identical whatever the number of threads. The reported "computation time" is wall-clock time spent in the simulation loop.

## Alternative

TROLL 4.0 can be set-up and run, and its outputs can be analyzed with an updated version of the R package rcontroll (https://github.com/sylvainschmitt/rcontroll/tree/TROLLV4), also available in R through the command:

```R
devtools::install_github("sylvainschmitt/rcontroll", ref = "TROLLV4").
```

The package also provide examples of input files and outputs.

## Help wanted?

Please use issues for any question regarding TROLL.

## Core group

* Jérôme Chave, CNRS, jerome.chave@univ-tlse3.fr
* Isabelle Maréchaux, INRAE, isabelle.marechaux@inrae.fr
* Fabian Fischer, TUM, f.j.fischer@tum.de
* Sylvain Schmitt, CIRAD, sylvain.schmitt@cirad.fr

