# Contrinbuting

TROLL is a scientific code where precisely identifying versions is crucial. Contributions follow the [Gitflow](https://nvie.com/posts/a-successful-git-branching-model/) model.

- `dev` is the default working branch.  
- `main` is protected (no direct commits). It only receives **merge commits** and **annotated tags** (`v3.9`, `v4.0`, `v4.0.1`, ...). Each minor release is documented in the `NEWS.md`.
- Everything starts with an **issue**, labeled as a `new feature`, `bug fix` or `hot fix`. Each issue gets a dedicated branch, named after the issue.
  - `new feature` and `bug fix` issues: branch from `dev`.
  - `hot fix` issues (unexpected bug in a released version): branch from `main`. 
- Branches are merged through **pull-requests** (from a fork or the same repository), which allows for review and CI/CD checks. Pull-requests targeting `dev` are merged when the issue is solved; those targeting `main` are merged when a release is ready.
- Long-term planning occurs in **milestones** named after upcoming minor versions (e.g. `3.9`, `4.0`).
- When a milestone is completed, `dev` is merged into `main` through a pull-request, and the minor version is tagged with an annotated tag (e.g. `v4.1`).
- Hot fixes are merged into `main` through a pull-request, tagged as a **patch version** (e.g. `v4.0.1`), **and merged back into `dev`**.
- Core team members work on branches directly in the TROLL repository, while external contributors work on a **fork** and open a pull-request from their fork. In both cases, pull-requests target `dev` (or `main` for hot fixes only).
