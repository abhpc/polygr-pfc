# Repository Guidelines

## Scope and scientific compatibility

Follow `.specify/memory/constitution.md`. Preserve the original physical equations,
discrete algorithms, FP64 arithmetic, input compatibility and two-stage graphene
workflow. Limit solver changes to DCU acceleration and necessary supporting work.
Do not rebuild the retained Coordinator JAR from newer sources with a different interface.

## Project layout

- `src/pfc.c`: shared solver driver; `src/pfc_dcu.hip.cpp`: HIP/hipFFT and MPI backend.
- `source/`: auxiliary C/Java tools, Python grayscale exporter and OVITO conversion entry point.
- `java/`: retained workflow JARs.
- `example/`: two-stage input templates.
- `slurm/`: dependency build, solver build, calculation scripts and PFC modulefile.
- `tests/`: minimal validation helper and frozen inputs.
- `specs/001-dcu-pfc-acceleration/`: specifications, contracts and anonymous validation results.

## Build and execution

Public users may compile and run directly with the shell commands in `README.md`;
Slurm and Debug are not prerequisites for using the software. Obtain compiler and
library paths from the installed toolchain.

For agents working in this managed session, all compilation, syntax checks, tests,
calculations and postprocessing must still run in the managed cluster's Slurm Debug
batch jobs. Never execute them locally or directly on its login node, submit to
another partition, or fall back. This agent/deployment restriction does not apply
to public users on their own systems.

The optional `slurm/fftw.slm`, `slurm/build.slm` and `slurm/polygr.slm` are
managed-deployment examples. For that deployment, register the repository module
directory and load `pfc/2026.1`; adapt module names to the site.

Keep `.slm` files as sequential commands with `#!/bin/bash -e` and resources in
`#SBATCH` directives. Do not add conditionals, workflow wrappers, `set -euo pipefail`,
`${VAR:?}` guards, repeated resource checks, outer `time`, probes, tests or evidence
collection to daily scripts. Preserve the solver's own per-step timing logs.

Both host builds use `-O3 -march=native -ffp-contract=off -Wall`; device arithmetic
is FP64 with contraction disabled. Link hipFFT and FFTW against the same MPI ABI.
The validated decomposition uses 32 total MPI ranks and four DCUs per node;
two-node runs distribute those ranks as 16 per node. Do not replace this with four
MPI ranks or bind every rank to an exclusive device. Preserve independent run
paths and save required project paths before changing module environments.

## Validation

Agents run only minimal validation in managed Debug batch jobs; no full suite, production-length
validation or parameter scan. Reuse existing results for documentation-only changes.
V1 uses the frozen short input; V2 reuses its CPU reference across two nodes.
V3 is the fixed 100×100 nm performance case: 4240² grid, 2000 steps per stage.
Measure complete two-stage PFC execution, excluding image and LAMMPS conversion.

Compare against original commit `c86efb51cf03086fb98897ab8abfe2704a62096e` with fixed
seed, MPI decomposition, physical parameters and output cadence. Keep the CPU source
unchanged. Field tolerances are `1e-6 + 1e-5*abs(CPU)`; scientific diagnostic
tolerances are `1e-8 + 1e-6*abs(CPU)`. Do not relax tolerances after seeing results.
A successful process or Slurm exit alone does not establish numerical agreement or speedup.
Keep detailed job provenance private and publish only anonymized validation summaries.

## Style and contributions

Preserve tab indentation and existing braces in C/Java. Use existing C `snake_case`,
Java `PascalCase` classes and `camelCase` methods. Keep formatting changes local;
compile C with `-Wall`. Use concise commit subjects. PRs should explain behavior,
validation, dependencies and numerical impact; include images for visualization changes.

## Public and private content

Keep public files free of account names, SSH endpoints, internal addresses, personal
absolute paths and private job records. Store site-specific connection details,
configuration copies and run artifacts in the Git-ignored `yhmsr/` directory.
Do not force-add that directory or generated outputs. Preserve upstream author
credits and license notices. Never replace an existing site's deployed configuration
with a generic template without adapting its paths and modules.
