# Implementation Plan: [FEATURE]

**Branch**: `[###-feature-name]` | **Date**: [DATE] | **Spec**: [link]

**Input**: Feature specification from `/specs/[###-feature-name]/spec.md`

**Note**: This template is filled in by the `/speckit-plan` command; its definition describes the execution workflow.

## Summary

[Extract from feature spec: primary requirement + technical approach from research]

## Technical Context

Public users may build and run with the direct shell commands in README; do not require
Slurm, Debug or a private cluster module for public installation. The execution restrictions
below apply to agents in this managed session and its private cluster deployment only.

<!--
  ACTION REQUIRED: Replace the content in this section with the technical details
  for the project. The structure here is presented in advisory capacity to guide
  the iteration process.
-->

**Language/Version**: [e.g., Python 3.11, Swift 5.9, Rust 1.75 or NEEDS CLARIFICATION]

**Primary Dependencies**: [e.g., FastAPI, UIKit, LLVM or NEEDS CLARIFICATION]

**Storage**: [if applicable, e.g., PostgreSQL, CoreData, files or N/A]

**Testing**: [CPU/DCU numerical equivalence, two-stage/conversion regression and performance;
agent builds and tests via managed Debug Slurm jobs; public shell use is described in README]

**Target Platform**: [Compatible MPI/DCU environment for public users; managed agents use
allocated Debug Slurm nodes and confirm devices, toolchain and FFT library]

**Project Type**: [e.g., library/cli/web-service/mobile-app/compiler/desktop-app or NEEDS CLARIFICATION]

**Performance Goals**: [domain-specific, e.g., 1000 req/s, 10k lines/sec, 60 fps or NEEDS CLARIFICATION]

**Constraints**: [domain-specific, e.g., <200ms p95, <100MB memory, offline-capable or NEEDS CLARIFICATION]

**Scale/Scope**: [domain-specific, e.g., 10k users, 1M LOC, 50 screens or NEEDS CLARIFICATION]

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

- [ ] **DCU scope**: Every code change directly enables acceleration or its required build,
  scheduling, compatibility or validation; unrelated work is excluded.
- [ ] **Physics unchanged**: Map migrated operations to the original CPU formulas and discrete
  algorithms; preserve FP64, dynamics, FFT normalization, boundaries and box optimization.
- [ ] **Compatibility**: Preserve inputs, CLI, restart/state data, both stages and LAMMPS conversion.
- [ ] **Managed-agent execution**: Agent compilation and tests, including CPU baseline and auxiliary
  validation, use sequential `.slm` scripts submitted through the target Slurm cluster with `sbatch path/to/job.slm`.
  Set `#SBATCH --partition=Debug` and resources in the script, use `#!/bin/bash -e` and explicit
  paths, and do not add `set -euo pipefail`, `${VAR:?}`, shell partition/node/rank guards
  or workflow wrappers.
  Never use another partition or fallback, local or login-node execution.
  Daily scripts contain only build/run commands: no tests, evidence/metadata generation,
  `/usr/bin/time`, probes or debug/validation calls. Preserve historical validation evidence
  separately and retain the solver's per-step timing logs.
- [ ] **Evidence**: Identify the CPU baseline, fixed inputs/initial fields/seeds, error metrics
  and justified tolerances before execution; plan numerical checks before matched-workload timing.
  Retain job IDs, resources, environment, exit status and logs.

These gates cannot be waived through Complexity Tracking. Unverified execution remains pending.

## Project Structure

### Documentation (this feature)

```text
specs/[###-feature]/
├── plan.md              # This file (/speckit-plan command output)
├── research.md          # Phase 0 output (/speckit-plan command)
├── data-model.md        # Phase 1 output (/speckit-plan command)
├── quickstart.md        # Phase 1 output (/speckit-plan command)
├── contracts/           # Phase 1 output (/speckit-plan command)
└── tasks.md             # Phase 2 output (/speckit-tasks command - NOT created by /speckit-plan)
```

### Source Code (repository root)
<!--
  ACTION REQUIRED: Replace the placeholder tree below with the concrete layout
  for this feature. Delete unused options and expand the chosen structure with
  real paths (e.g., apps/web, packages/something). The delivered plan must
  not include Option labels.
-->

```text
# [REMOVE IF UNUSED] Option 1: Single project (DEFAULT)
src/
├── models/
├── services/
├── cli/
└── lib/

tests/
├── contract/
├── integration/
└── unit/

# [REMOVE IF UNUSED] Option 2: Web application (when "frontend" + "backend" detected)
backend/
├── src/
│   ├── models/
│   ├── services/
│   └── api/
└── tests/

frontend/
├── src/
│   ├── components/
│   ├── pages/
│   └── services/
└── tests/

# [REMOVE IF UNUSED] Option 3: Mobile + API (when "iOS/Android" detected)
api/
└── [same as backend above]

ios/ or android/
└── [platform-specific structure: feature modules, UI flows, platform tests]
```

**Structure Decision**: [Document the selected structure and reference the real
directories captured above]

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| [e.g., 4th project] | [current need] | [why 3 projects insufficient] |
| [e.g., Repository pattern] | [specific problem] | [why direct DB access insufficient] |
