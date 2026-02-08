# CLAUDE.md

## Project Overview

**LM** is a numerical solver for Ornstein-Zernike (OZ) equations with Percus-Yevick (PY) and Hypernetted Chain (HNC) closures in two dimensions. It handles both monodisperse (single-component) and binary mixture systems. Written in C99, originally developed in 2005 as part of a diploma thesis in statistical mechanics.

**Author**: Gerolf Ziegenhain
**License**: GPL v3

## Repository Structure

```
lm/
├── CLAUDE.md              # This file
├── COPYING                # Copyright notice
├── README.txt             # Brief project description
├── TechRep.pdf            # Technical report (mathematical background)
├── gpl-3.0.txt            # GPL v3 license text
├── monodisperse/          # Single-component solver
│   ├── makefile
│   ├── static.c           # Main program (~1,700 lines)
│   └── marquardt.c        # Marquardt fitting algorithm (~1,400 lines)
└── binary/                # Binary mixture solver
    ├── makefile
    ├── static.c           # Main program (~2,300 lines)
    ├── marquardt.c        # Marquardt fitting algorithm (~1,400 lines)
    └── PlotCorrelators.sh # Gnuplot visualization script
```

The two subdirectories (`monodisperse/` and `binary/`) are independent programs with near-identical `marquardt.c` fitting code (see differences below).

## Build System

Each subdirectory has its own `makefile`. Build with:

```sh
cd monodisperse && make    # builds ./static
cd binary && make          # builds ./static
make clean                 # removes built executable
```

**Compiler flags**: `-Wall -pedantic -funroll-loops -std=c99`
**Link flags**: `-lgsl -lgslcblas -lm -O3`

The output executable is named `static` in each directory.

### Dependencies

| Dependency | Purpose |
|---|---|
| **GCC** | C99 compiler |
| **GSL** (GNU Scientific Library) | Bessel function zeros (`gsl_sf_bessel_zero_J0`) |
| **GSL CBLAS** | Linear algebra (linked via `-lgslcblas`) |
| **libm** | Standard C math |
| **lmfit** | Levenberg-Marquardt fitting (bundled in `marquardt.c`) |

GSL is optional — there is a fallback formula for Bessel zeros when GSL is unavailable.

## Code Architecture

### Computational Pipeline

1. **Potential calculation** (`return_mayer()`) — computes Mayer f-functions for the interaction potential
2. **Closure relation** — PY or HNC approximation (selected via `#define APPROXIMATION`)
3. **OZ equation solver** (`OZ()`) — core integral equation solver
4. **Fast Bessel Transform** (`FBT1()` / `FBT2()`) — converts between real-space and reciprocal-space
5. **Marquardt optimizer** (`lm_minimize()`) — non-linear least squares fitting
6. **Interpolation** — 2D spline interpolation (`spline()` / `splint()`)

### Key Data Structures

```c
typedef struct { unsigned int dim; double **data; } t_matrix;
typedef struct { unsigned int step; t_matrix data; } t_fbt;
```

### Configuration via Preprocessor Defines

Parameters are set at compile time in `static.c`:

| Define | Values | Meaning |
|---|---|---|
| `SYSTEM` | 1–3 | Interaction potential type (e.g., 3 = dipolar) |
| `APPROXIMATION` | 1 or 2 | 1 = PY, 2 = HNC |
| `STARTVALUE` | 1 or 2 | Initial condition for Gamma(r) |
| `BASIS_LINE` | defined/undefined | Use line-based basis decomposition |
| `PACKING_FRACTION` | defined/undefined | Specify packing fraction eta |
| `SYMMETRIC_FT` | defined/undefined | Symmetric Fourier transform |
| `MARQUARDT` | defined/undefined | Enable Marquardt optimizer |

### Debug/Test Modes

Compile-time testing is controlled via defines:

- `#define DEBUG 100` — limits iterations for quick debugging
- `#define TEST N` — runs specific test routines (N = 1–7):
  - 1: Fast Bessel Transform, 2: Decomposition, 3: Start Gamma, 5: Inversion, 6: Real decomposition, 7: Gamma reading

### Output Files

The solver writes correlation function data to text files:

- `Cr__pic*.dat` — direct correlation function C(r)
- `Hr__pic*.dat` — total correlation function H(r)
- `Gammar__pic*.dat` — bridge function Gamma(r)
- `Cq__pic*.dat`, `Hq__pic*.dat` — reciprocal-space equivalents

## Development Notes

- **No formal test suite** — testing is done via compile-time `#define TEST` modes
- **No CI/CD** — this is legacy academic code
- **No linter config** — relies on `-Wall -pedantic` compiler warnings
- **Comments are partially in German** — legacy from the original diploma thesis
- **No package manager** — dependencies are system-level (install GSL via your OS package manager)
- The `marquardt.c` files are **not identical** between `monodisperse/` and `binary/`. Key differences:
  - `lm_evaluate_default` declaration: commented out in monodisperse, active in binary
  - `#ifdef DEBUG` output: commented out in monodisperse, active in binary
  - **Convergence tolerance**: `epsilon = 1.e-14` (monodisperse) vs `1.e-10` (binary) — factor 10,000 difference in Marquardt precision

## Conventions for AI Assistants

- This is pure C99 code. Do not introduce C++ features or non-standard extensions.
- Preserve the existing code style: procedural C with minimal abstraction.
- Configuration changes are made by editing `#define` values in `static.c`, not via command-line arguments.
- The two solver directories are independent — changes to shared logic (e.g., `marquardt.c`) should be applied to both copies.
- Refer to `TechRep.pdf` for the mathematical background behind the algorithms.
- When building, ensure GSL development headers are installed (`libgsl-dev` or equivalent).
