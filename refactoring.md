# Refactoring notes

> This documents the refactoring process of `src/troll.cpp` from a ~10 000-line monolith into a modular, context-threaded codebase.

---

## Starting point: the monolith

Before this refactor, the entire simulation lived in a single file. It had ~10 000 lines and ~200 file-scope global variables — parameters, field arrays, output streams, RNG state, and the `Species`/`Tree` vectors were all declared at file scope.

---

## The quick dev-test: `make sanity`

Before any code was moved, a regression test was set up: `make sanity` compiles both the original TROLL v4.0 monolith and the refactored version, runs both on the same inputs, and compares the 41 output files. 

Lookup tables were cached for repeated checks to speed things up.

This quick-check (20 iterations) was run after every commit to ensure nothing "broke". 

In the end 1000 iterations were tried too successfully. 

---

## Phase 1 — Constant/flag extraction (the easy part)

Initially, three small headers were split out:
- `include/troll_defines.hpp` — compile-time feature flags (`#define`/`#undef WATER`, `Output_ABC`, `CACHE_LUT`, etc.)
- `include/constants.hpp` — physical constants (`PARSEC`, `kBoltzmann`, etc.)
- `include/troll.hpp` — the remaining global declarations (later reduced further)

---

## Phase 2 — The `Context` struct

A single `Context` struct was introduced in `include/context.hpp`, grouping all ~200 globals into 17 named sub-structs:

| Sub-struct | Contents |
|---|---|
| `FileIO` | Input/output file path char buffers |
| `InputBuffers` | Raw `argv` pointers for command-line parsing |
| `OutputConfig` | All `fstream` output file objects |
| `ModelOptions` | Boolean feature flags (`_BASICTREEFALL`, `_NDD`, etc.) |
| `Grid` | Site count, dimensions, spatial arrays, HEIGHT, SBORD |
| `TimeState` | `iter`, `nbiter`, `iterperyear`, `timestep`, etc. |
| `Climate` | Daily climate vectors, per-step scalars, yearly means |
| `LookupTables` | All `LookUp_*` pointer arrays and their bin counts |
| `Soil` | ~25 soil/water variables and 3-D field arrays |
| `Diagnostics` | Live/dead tree counts, DBH histograms |
| `PointCloud` | LiDAR simulation parameters |
| `RNGState` | GSL RNG pointer, covariance matrix/vectors |
| `SpeciesState` | Seed arrays, species probabilities |
| `Intraspecific` | 10×10000 trait deviation arrays, sigmas, correlations |
| `CrownGeometry` | Gap fraction, crown shape, visual extent bounds |
| `SimParams` | 93 simulation parameters (light, CO2, allometry, mortality…) |
| `SimFields` | `LAI3D`, `Thurt[3]` dynamic field arrays |

`vector<Species> S` and `vector<Tree> T` are direct fields on `Context`, not wrapped in a sub-struct (see Phase 5).

**Grouping?**

Being unaware of the scientific context, quite some guessing took place here: variables that are always read and written together, or that form a concept, belong together. Hope it makes sense...

The grouping also controls scope. A function that only queries the grid does not need to know that `Soil` or `RNGState` exist.

**A note on performance**

`OutputConfig` is placed **last** in `Context`. `fstream` objects are big and rarely accessed during the simulation. Placing them last keeps the frequently-accessed fields (`grid`, `time`, `params`, `field`) at low offsets for better caching, otherwise it was observably slower.

---

## Phase 3 — Splitting `troll.cpp`

`troll.cpp` was split into eight translation units:

| File | What it contains | Lines |
|---|---|---|
| `src/troll.cpp` | `main()` only — orchestration | ~318 |
| `src/species.cpp` | `Species::Species()`, `Species::Init()` | ~54 |
| `src/tree.cpp` | All `Tree::*` methods | ~3355 |
| `src/crown.cpp` | Crown geometry and light calculation helpers | ~353 |
| `src/lookup_tables.cpp` | `InitialiseLookUpTables`, LAImax, intraspecific LUTs | ~414 |
| `src/input.cpp` | All `ReadInput*` functions, species/pointcloud assignment | ~842 |
| `src/memory.cpp` | `AllocMem`, `FreeMem`, `Initialise`, output stream init, ABC init | ~1183 |
| `src/output.cpp` | `Average`, all `Output*` functions, close functions | ~3189 |
| `src/simulation.cpp` | `Evolution`, `UpdateSeeds`, `UpdateField`, `FillSeed`, `RecruitTree`, treefalls | ~745 |

**Grouping?**

Again, fingers crossed for my splitting decisions.

---

## Phase 4 — Threading `Context` by reference

`Context ctx` was made a local variable in `main()`. Every function that previously touched a global was updated to accept `Context &ctx` as its first argument and access state through it.


This threading was done, function-group by function-group (`crown`, `lookup`, `species`, then `input`, then `memory`/`output`/`simulation`). Each was independently sanity-tested.

---

## Phase 5 — Species and Tree arrays

- `Species` class declaration was extracted from `troll.hpp` into `include/species.hpp`.
- `Tree` class declaration (and the `leafFluxes` struct) was extracted into `include/tree.hpp`.
- `include/context.hpp` was updated to include: `vector<Species> S` and `vector<Tree> T` as direct `Context` fields.
- All `S[i]` / `T[i]` call sites across all `.cpp` files were replaced with `ctx.S[i]` / `ctx.T[i]`.


---

## The parameter registry

The 300-line `AssignValueGlobal` if/else chain (which assigned values to globals after reading `global_inputs.txt`) was replaced by a registry in `src/params/param_registry.cpp`. 

Each parameter is registered once with its name, target pointer, min value, max value, and default value. 

A single `SetParameter<T>` template validates and assigns. 

The same pattern was later applied to species parameters (`BuildSpeciesRegistry` / `AssignSpeciesParam` in `src/input.cpp`).

**A note: Context before input**

The registry binds directly to `ctx.*` fields via pointers, so it requires a `Context &ctx` at registration time. `RegisterParameters(ctx)` must be called early in `main()`, before any file is read.

## Summary

| Dimension | Before | After |
|---|---|---|
| Largest single file | ~10 000 lines (`troll.cpp`) | ~3355 lines (`tree.cpp`) |
| Global variable count | ~200 file-scope | 3 (MPI, disabled by default) |
| Parameter assignment | 300-line if/else chain | Registry, 93 entries |
| Function signatures | No explicit state | `Context &ctx` first argument |
| Multi-TU compilation | Single translation unit | 9 translation units, parallel build |
| Testability | One global simulation state | `Context` is a local value; two contexts = two simulations |
| Data flow | Implicit (any function, any global) | Explicit (visible in call site) |

---

## TODOs?

### MPI globals

`mpi_rank`, `mpi_size`, and `easympi_rank` remain as file-scope globals in `troll.cpp`. MPI support is disabled by default (`#undef MPI`). Moving these into `Context` would require threading them through functions that are not currently MPI-aware.

In a conversation here [https://github.com/troll-model/TROLL/pull/22#issuecomment-3580696619] I saw that maybe it should be removed? TBD

### Conditional-compilation diagnostics

`#ifdef Output_ABC` and `#ifdef CHECK_CARBON` blocks remain conditional. The ABC functions were updated to accept `Context &ctx`, but the ABC globals (`chm_field_*`, `transmittance_*`, moving averages) are **not** in `Context` yet

### Local variable renames (marked  with `// RENAMED:` comments)

During the split (Phase 3) and the `Context` threading (Phase 4), several local variables and function parameters were renamed to avoid shadowing same-named `ctx.params.*` fields. Each renamed site carries a `// RENAMED:` comment explaining the original name and why it changed.

The affected names are:

| Original name | New name | Location | Reason |
|---|---|---|---|
| `alpha` | `alpha_vgm` | `src/input.cpp` | Van Genuchten–Mualem soil parameter; distinct from `ctx.params.alpha` (apparent quantum yield) |
| `alpha` | `alpha_crown` | `src/tree.cpp` | FvCB canopy integration ratio; distinct from `ctx.params.alpha` |
| `dens` | `dens_layer` | `src/crown.cpp` (params + locals) | Per-layer crown density; was shadowing `ctx.params.dens` |
| `dens` | `dens_lowerlayer` | `src/crown.cpp` (locals) | Base-layer density from `GetDensitiesGradient`; distinct from the `dens_layer` output parameter |

**We can probably revert back if needed** 

After threading `Context`, the former globals are accessed as e.g. `ctx.params.alpha`, a local `alpha` can no longer shadow them. 

---

### Confidence ###

### The scientist's validation ###

So far the only test that I haven't broken anything is the sanity check described above. It does compare row by row, so per-step divergencies are checked, as well as final values. I tried it at 20 iterations while developing and at 1000 once when I finished. Perhaps a bigger number should be tried? Not sure how the model is behaving in time. Or with different inputs/initial params etc. 

Most certainly, you, if you are reading this, will have more ideas than me on how to better test the "science-correctness" of this refactoring.


### The developer's validation ###

Compilation warnings are not worrying (one could fix/remove them though). So, in combination with the sanity checks performed and the knowledge that I mostly "moved things around" I am quite confident this is sound. 


### A note on a tricky bug (hope there's no more of these that just didn't come up) ###

During the refactoring, one bug was very tricky to find: 
CloseOutputs only looped i < 10, leaving output[11], output[31], output[32] (and others) never closed. Those were flushed automatically when `exit(0)` triggered global destructor cleanup — because `ctx` was global. With `ctx` local, `exit(0)` skiped local destructors, so those buffers were silently lost. Perhaps we should check if similar things exist elsewhere.


