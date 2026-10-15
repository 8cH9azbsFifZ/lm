## Why

In October 2026 the diploma thesis (2006) was recomputed from scratch. The work had three parts:

- The original code was rebuilt and rerun. It reproduces the thesis data to 1e-10, and binary hard discs bit-identically.
- The results were checked independently against a separate Ornstein–Zernike solver, an event-chain Monte Carlo simulation, and two independent mode-coupling implementations.
- The code was reviewed.

The static correlators are confirmed. The review found defects that can produce wrong results or wrong reports:

1. **False convergence (monodisperse).** The convergence measure is divided by the signed sum of the fine part. When that sum is negative, the measure is negative, and the loop reports convergence after a few iterations. In 2026 this happened in 3 of 4 strongly coupled dipolar test runs.
2. **Type punning (monodisperse).** A `long double` array is passed to the Levenberg–Marquardt routine as `double*`. On x86 (80-bit long double) the LM start values are bit garbage.
3. **Out-of-bounds read (both solvers).** The last basis function reads `i_[nu]`, one element past the end of the array.
4. **Wrong Γ diagnostic (monodisperse).** The printed Γ uses `10e-7` instead of `1e-7` (μ0/4π), a factor of 10. It also passes `long double` values to `%g`.
5. **Binary solver does not compile for SYSTEM 2.** `delta` is redefined under `PACKING_FRACTION`.
6. **Dipolar parameters (monodisperse) do not match the thesis data.** The thesis data were computed with σ = 2.8 µm (r = 1.4 µm) and χ = 6.2e-12. The code had r = 2.35 µm.
7. **MCT quadrature bias.** In the 2D vertex, the integrable 1/√ singularities at p = |q−k| and p = q+k are sampled at mid-cell points with Simpson end weights, and the last cell is dropped. The vertex is underestimated, and η_c^PY comes out as 0.7044 instead of 0.6915 (literature: 0.6914, Hajnal, Brader & Schilling, PRE 80, 021503 (2009)). The MCT code was not in this repository until now.

## What Changes

- Normalise the monodisperse convergence measure by |Σ fine|. The result is bit-identical whenever Σ fine > 0, which covers all documented 2006 runs.
- Pass a separate `double` buffer to `lm_minimize` and copy the values back afterwards.
- Guard the last basis function against the out-of-bounds index.
- Correct the Γ diagnostic (1e-7, casts for printf).
- Fix the `delta` redefinition in `binary/static.c`.
- Set the monodisperse particle radius to r = 1.4 µm, as used for the thesis data, and document χ.
- Add an optional environment variable `LM_PICMAX` to cap the iterations (for tests).
- Add the MCT solver `mct/nep.c` (2006 listing) with a compile switch `QUAD`:
  - `QUAD=1` (default): exact cell integration of the singular factor.
  - `QUAD=0`: the original 2006 scheme, for reproducing the thesis.
- Add `ERRATUM.md` and tests.

## Impact

- `monodisperse/static.c`, `binary/static.c`, new `mct/`, new `tests/`, `README.txt`, new `ERRATUM.md`.
- No change to converged static results. The thesis value η_c = 0.7044 is superseded by 0.6915.
