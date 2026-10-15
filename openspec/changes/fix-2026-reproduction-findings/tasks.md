## 1. Static solver (monodisperse)
- [x] 1.1 Normalise the convergence measure by |Σ fine|
- [x] 1.2 Add a `double` buffer for `lm_minimize`
- [x] 1.3 Guard the last basis function in both basis paths against `i_[nu]`
- [x] 1.4 Correct the Γ diagnostic (1e-7, casts)
- [x] 1.5 Set r = 1.4e-6 and document χ
- [x] 1.6 Add the optional iteration cap `LM_PICMAX`

## 2. Static solver (binary)
- [x] 2.1 Fix the `delta` redefinition under `PACKING_FRACTION`
- [x] 2.2 Guard the last basis function against `i_[nu]`

## 3. MCT
- [x] 3.1 Add `mct/nep.c` (2006 listing) with `QUAD` switch and start value f = 1
- [x] 3.2 Add the makefile

## 4. Tests and documentation
- [x] 4.1 Add fixtures from the thesis and from the independent PY solver (`tests/data/`)
- [x] 4.2 Add `tests/run_tests.py` covering the spec scenarios
- [x] 4.3 Add `ERRATUM.md` and update `README.txt`
- [x] 4.4 Verify on the target system (V1–V5) and record the results in ERRATUM.md
