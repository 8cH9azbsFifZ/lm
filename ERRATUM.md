# Erratum and recomputation (October 2026)

The diploma thesis (2006) and this code were recomputed from scratch. The work had three parts:
- The original code was rebuilt on current hardware (Apple Silicon, clang, GSL 2.8) and run again.
- The results were checked independently against:
  - a separate 2D Ornstein–Zernike solver (Newton–Krylov on a uniform grid),
  - an event-chain Monte Carlo simulation of hard discs (N = 4096 and 16384),
  - two independent mode-coupling implementations.
- The code was reviewed line by line.

## What holds

- **Static correlators are reproducible.**
  - Hard discs, monodisperse: n σ² = 0.71, PY, matches the thesis data to 2·10⁻¹⁰.
  - Hard discs, binary: η = 0.5, δ = 0.5, bit-identical (3·10⁻¹⁴), with the same iteration count as the 2006 cluster run.
- **The independent OZ solver confirms the correlators.** Residual differences of 10⁻³ to 10⁻² come from the Bessel-zero grid: there is no node at r = σ, so the effective diameter is 0.997–1.004 σ.
- **The claim that Picard iteration fails in 2D at high density holds.** Above η ≈ 0.67 nothing converges, so the Levenberg–Marquardt/Gillan approach is justified.

## Corrections to the thesis

| Thesis | Correct | Cause |
|---|---|---|
| η_c^PY = 0.7044 (MCT, monodisperse hard discs) | **0.6915 ± 0.0005** | Quadrature bias in the 2D vertex (see below). Literature: 0.6914, Hajnal, Brader & Schilling, PRE 80, 021503 (2009) |
| η_c^MHNC = 0.712 | ≈ 0.700 (estimate) | Same bias. Literature: 0.697, Bayer et al., PRE 76, 011508 (2007) |
| Double peak in the dipolar NEP is "independent of discretisation" | At η = 0.67 only a shoulder remains | Same bias |
| PY data at n σ² ≥ 0.962 | Not converged (stopped at 100000 iterations); S(q) peak about 15 % too low | Iteration limit |
| PY becomes worse at high density "because it is a linearisation of HNC" | PY is clearly better than HNC at all densities (MC: contact value PY −2 … −11 %, HNC +3 … +49 %) | Interpretation |
| Dipolar parameters: χ = 6.6 pAm²/T | Data were computed with χ = 6.2·10⁻¹², r = 1.4 µm (βV(σ)·T = 12683.1 K) | Documentation |
| Γ_m values in Ch. 3/6, e.g. Γ_m^c = 2…5 | About 4× larger by definition (eq. gamma_mono). The transition line in (η, T) is unchanged | Normalisation unclear |

**Context (not a mistake in the thesis):** For monodisperse hard discs, η ≈ 0.70–0.72 lies in the liquid–hexatic coexistence region (Bernard & Krauth, PRL 107, 155704 (2011)).

### MCT quadrature

The p integral in the 2D vertex contains the factor 1/√(k² − ((q²+k²−p²)/2q)²). This factor diverges, but integrably, at p = |q−k| and p = q+k.

`nep.c` (2006):
- sampled this factor at cell midpoints with Simpson end weights,
- dropped the last cell next to q+k.

As a result the vertex is systematically too small, with a bias proportional to √dq. On identical, accurate PY structure factors:

| η | 0.691 | 0.692 | 0.704 | 0.706 |
|---|---|---|---|---|
| `QUAD=1` (exact cell integration) | liquid | **glass** | glass | glass |
| `QUAD=0` (2006 scheme) | liquid | liquid | liquid | **glass** |

## Code defects fixed in 2026

1. **False convergence (monodisperse).** The convergence measure was divided by the *signed* sum of the fine part. When that sum was negative, the run reported "converged".
   - In 2026 this occurred with strong dipolar coupling.
   - It does not occur in any of the 32 preserved 2006 logs.
   - Now divided by |Σ|. This is bit-identical to before whenever Σ > 0.
2. **Type punning.** `lm_minimize` received a `long double*` as `double*`. On x86 the start values were bit garbage. A `double` buffer is now passed instead. The fixed point is unchanged.
3. **Out-of-bounds read.** `i_[nu]` was read for the last basis function (both solvers).
4. **Γ diagnostic** used `10e-7` instead of `1e-7`. This affected the screen output only.
5. **`binary/static.c` did not compile for SYSTEM 2** (`delta` was redefined).
6. **Monodisperse radius.** It was r = 2.35 µm; the thesis data need r = 1.4 µm.
7. **New `mct/` directory.** The MCT solver was added with the switch `QUAD` (1 = corrected, 0 = 2006 scheme).

Note: With the thesis coupling (e.g. n σ² = 0.5, T = 1300 K, B = 1 mT), the static solver does not converge from the start value γ = 0. Such cases need a start value from a weaker coupling (`STARTVALUE 3`).

## Tests

```sh
python3 tests/run_tests.py          # fast
python3 tests/run_tests.py --full   # all scenarios, about 45 min
```

The scenarios are specified in `openspec/changes/fix-2026-reproduction-findings/`.

## Verification on the target system (2026-10-15)

`python3 tests/run_tests.py --full` (Apple M-series, macOS, GSL 2.8): **10/10 passed**

| Scenario | Result |
|---|---|
| V1 HS-PY n σ² = 0.71 vs thesis S(q) | max\|ΔS\| = 1.2·10⁻¹⁰ |
| V2 no negative convergence measure, unconverged → −1 | OK |
| Calibration βV(σ)·T | 12683.103 K |
| V3 binary compiles for SYSTEM 2/3 | OK |
| V4 `QUAD=0` vs thesis fq_0.9 | max\|Δf\| = 1.3·10⁻⁹ |
| V5 `QUAD=1` η = 0.691 / 0.692 | liquid / glass (f_max = 0.867) |
