## ADDED Requirements

### Requirement: Non-ergodicity parameters with exact singular integration
The MCT solver `mct/nep.c` SHALL integrate the integrable 1/√ singularities of the 2D vertex exactly over each cell (`QUAD=1`, the default).

#### Scenario: Liquid below the critical density
- **WHEN** `nep` runs with QUAD=1, 500 points and qmax = 50 on the PY S(q) at η = 0.691
- **THEN** the iteration ends at f ≡ 0 (norm < 1)

#### Scenario: Glass above the critical density
- **WHEN** `nep` runs with QUAD=1, 500 points and qmax = 50 on the PY S(q) at η = 0.692
- **THEN** the iteration ends at a nonzero f with norm > 1

### Requirement: Reproducing the thesis scheme
With `QUAD=0` the MCT solver SHALL reproduce the 2006 discretisation.

#### Scenario: Thesis NEP file
- **WHEN** `nep` runs with QUAD=0, 500 points, qmax = 50 and the thesis S(q) `sq_0.9`
- **THEN** f(q) deviates from the thesis file `fq_0.9` by at most 1e-6
