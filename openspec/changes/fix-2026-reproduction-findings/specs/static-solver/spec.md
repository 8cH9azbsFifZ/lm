## ADDED Requirements

### Requirement: Converged correlators unchanged
For every converged case, the static solver SHALL deliver the same correlators as the 2006 version.

#### Scenario: Monodisperse hard discs reproduce the thesis data
- **WHEN** `monodisperse` is run with SYSTEM 1, PY, n σ² = 0.71, Rmax = 20, M = 1024 and the 2006 basis
- **THEN** S(q) deviates from the thesis file `hs_mono/py/0.71/Sq.dat` by at most 1e-8

### Requirement: No false convergence
The monodisperse solver SHALL report convergence only when the convergence measure is non-negative and below the threshold.

#### Scenario: Negative fine-part sum
- **WHEN** the sum of the fine part becomes negative during the iteration, as in the dipolar case n σ² = 0.5, T = 1300 K, B = 1 mT, M = 500, Rmax = 25
- **THEN** no negative `fine:` value is printed
- **AND** the run does not end as converged, as long as the measure is above the threshold

#### Scenario: Iteration cap for tests
- **WHEN** the environment variable `LM_PICMAX` is set to N
- **THEN** the solver stops after at most N iterations and returns -1 if it has not converged

### Requirement: Type-correct LM interface
The solver SHALL pass `double` values to `lm_minimize`, so that the start values are correct independent of the size of `long double`.

#### Scenario: LM start values on x86
- **WHEN** the code is built on a platform where `sizeof(long double) > sizeof(double)`
- **THEN** `lm_minimize` receives the current coarse coefficients as `double`
- **AND** the result is copied back into the `long double` array afterwards

### Requirement: Memory-safe basis construction
The basis construction SHALL NOT access `i_[nu]`.

#### Scenario: Last basis function
- **WHEN** the basis function with index ν−1 is built
- **THEN** only its ascending part is set, with no access past the end of the array

### Requirement: Correct dipolar parameters and diagnostics
The monodisperse dipolar solver SHALL use the parameters of the thesis data and print Γ correctly.

#### Scenario: Calibrating the Mayer function
- **WHEN** SYSTEM 2 is run with T = 200000 K and B = 1 mT
- **THEN** the Mayer function satisfies −ln(1+f(r))·r³ = 0.0634155 (± 1e-6)

#### Scenario: Γ output
- **WHEN** the solver prints Γ
- **THEN** the value is μ0/4π · B² χ² (πρ)^{3/2} / (k_B T), with μ0/4π = 1e-7

### Requirement: Binary solver compiles for all documented systems
`binary/static.c` SHALL compile for SYSTEM 2 and SYSTEM 3.

#### Scenario: Build for SYSTEM 2
- **WHEN** `binary/static.c` is compiled with `-DSYSTEM=2`
- **THEN** compilation succeeds without errors
