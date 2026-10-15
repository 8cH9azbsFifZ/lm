## Context

The code is written in C99 and the build uses GSL. The goal is maximum backward compatibility: every documented 2006 run must give identical results, and only the defective paths change. Success has to show on real data. Where possible, the tests compare against the original thesis data (fixtures in `tests/data/`).

## Decisions

- **Convergence measure:** divide by `fabsl(norm_sum)` instead of `norm_sum`. Alternatives considered:
  - Σ|fine| changes the stopping point of every run, so it is not bit-compatible.
  - The maximum norm has the same problem.
  - `fabsl` changes nothing as long as Σ fine > 0. In that case the old measure was already correct.
- **LM buffer:** keep a `double coarse_d[nu]` in the main loop and copy before and after `lm_minimize`. On platforms with `long double == double` (ARM64) this is a no-op.
- **Last basis function:** only evaluate the descending ramp if `n+1 < nu`.
- **MCT quadrature:** with w = q²+k²−p² and A = 2qk:

  ∫_a^b p dp / √(k² − ((q²+k²−p²)/2q)²) = q [asin(w_a/A) − asin(w_b/A)]

  The rest of the integrand, S_p f_p V², is evaluated at the cell midpoint. The p cells j = |l−i| … l+i cover exactly [|q−k|, q+k]. On the midpoint grid the k sum uses the plain midpoint rule, without Simpson end weights.
- **Dipolar parameters:** use r = 1.4e-6 m and χ = 6.2e-12. These values reproduce the stored Mayer functions of the thesis exactly: βV(σ)·T = 12683.1 K. The thesis text says χ = 6.6e-12. That discrepancy is documented in ERRATUM.md, not "fixed", because the data were computed with 6.2.

## Verification on the target system (success criteria)

- **V1:** The monodisperse HS-PY case (n σ² = 0.71, M = 1024, Rmax = 20, 2006 basis) reproduces `hs_mono/py/0.71/Sq.dat` from the thesis to ≤ 1e-8.
- **V2:** The dipolar test case (n σ² = 0.5, T = 1300 K, B = 1 mT, M = 500, Rmax = 25) no longer reports convergence with a negative measure, and no negative `fine:` value appears.
- **V3:** `binary/` compiles for SYSTEM 2 and SYSTEM 3.
- **V4:** `mct/nep` with `QUAD=0` reproduces the thesis file `fq_0.9` to ≤ 1e-6.
- **V5:** `mct/nep` with `QUAD=1` gives a liquid state (f ≡ 0) at η = 0.691 and a glass at η = 0.692, using PY S(q) from an independent solver (fixtures).

## Risks

- On x86 the LM path changes because the start values are now correct. The fixed point is the same; only the number of iterations can differ.
