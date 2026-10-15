#!/usr/bin/env python3
"""Tests for the scenarios in openspec/changes/fix-2026-reproduction-findings.

Usage:
    python3 tests/run_tests.py           # fast scenarios (~1 min)
    python3 tests/run_tests.py --full    # all scenarios, including V1, V4 and V5 (~45 min)

Requirements: gcc, GSL (on macOS via Homebrew), Python 3 standard library only.
On macOS, j0l/j1l are missing from libm; a shim based on j0/j1 is compiled in.
"""
import glob
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "tests", "data")
BASIS_HS = ".25 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.5 4 4.5 5 6 7 8 9 15".split()
BASIS_DIP = (".2 .4 .65 .9 1.5 1.7 2.1 2.3 2.7 3.0 3.5 4.0 4.5 5.0 7.0 8.0 9.0 10.0 11.0 "
             "12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0").split()
SHIM = "#include <math.h>\nlong double j0l(long double x){return j0((double)x);}\n" \
       "long double j1l(long double x){return j1((double)x);}\n"


def cflags():
    f = ["-std=c99", "-O3", "-w", "-DM_PIl=3.141592653589793238462643383279502884L"]
    for p in ("/opt/homebrew", "/usr/local"):
        if os.path.isdir(os.path.join(p, "include", "gsl")):
            f += ["-I" + p + "/include", "-L" + p + "/lib"]
    return f


def build_static(work, subdir, defines, name):
    """Copy static.c, set the #define switches and compile."""
    src = open(os.path.join(ROOT, subdir, "static.c"), encoding="latin1").read()
    for key, val in defines.items():
        src, n = re.subn(r"^#define %s \d+" % key, "#define %s %d" % (key, val), src, count=1, flags=re.M)
        assert n == 1, key
    d = os.path.join(work, name)
    os.makedirs(d)
    open(os.path.join(d, "static.c"), "w", encoding="latin1").write(src)
    shutil.copy(os.path.join(ROOT, subdir, "marquardt.c"), d)
    files = [os.path.join(d, "static.c")]
    if platform.system() == "Darwin":
        open(os.path.join(d, "shim.c"), "w").write(SHIM)
        files.append(os.path.join(d, "shim.c"))
    exe = os.path.join(d, "static")
    subprocess.run(["gcc"] + cflags() + files + ["-lgsl", "-lgslcblas", "-lm", "-o", exe], check=True)
    return exe


def build_nep(work, quad):
    exe = os.path.join(work, "nep_q%d" % quad)
    subprocess.run(["gcc", "-std=c99", "-O3", "-w", "-DQUAD=%d" % quad,
                    os.path.join(ROOT, "mct", "nep.c"), "-lm", "-o", exe], check=True)
    return exe


def read2(path):
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            rows.append((float(p[0]), float(p[1])))
    return rows


def latest(pattern, cwd):
    files = glob.glob(os.path.join(cwd, pattern))
    return max(files, key=lambda s: int(re.findall(r"__(\d+)\.dat", s)[0]))


def run(cmd, cwd, env=None):
    e = dict(os.environ)
    e.update(env or {})
    p = subprocess.run(cmd, cwd=cwd, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode("latin1")


RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, ok))
    print("%-62s %s  %s" % (name, "OK  " if ok else "FAIL", detail), flush=True)


def t_binary_compiles(work):
    for s in (2, 3):
        try:
            build_static(work, "binary", {"SYSTEM": s, "APPROXIMATION": 1}, "bin_s%d" % s)
            check("V3 binary compiles for SYSTEM %d" % s, True)
        except subprocess.CalledProcessError:
            check("V3 binary compiles for SYSTEM %d" % s, False)


def t_dipolar_calibration_and_no_false_convergence(work):
    exe = build_static(work, "monodisperse", {"SYSTEM": 2, "APPROXIMATION": 1}, "mono_dip")
    cwd = os.path.join(work, "v2")
    os.makedirs(cwd)
    rc, out = run([exe, "0.5", "1300", "1", "25", "500"] + BASIS_DIP, cwd, {"LM_PICMAX": "40"})
    fines = [float(x) for x in re.findall(r"fine: (\S+)", out)]
    check("V2 no negative convergence measure", all(f >= 0 for f in fines), "min=%g" % min(fines))
    check("V2 unconverged run ends with -1", rc == 255, "rc=%d, %d iterations" % (rc, len(fines)))
    # Calibration: -ln(1+f) r^3 * T = 12683.1 K (sigma = 2.8 um, chi = 6.2e-12, B = 1 mT)
    r, f = next((r, f) for r, f in read2(os.path.join(cwd, "mayer__0.dat")) if r > 1.0001)
    gt = -math.log(1.0 + f) * r ** 3 * 1300.0
    check("Calibration of the dipolar Mayer function", abs(gt - 12683.103) < 0.01, "betaV(sigma)*T=%.3f K" % gt)
    g = float(re.search(r"Gamma: (\S+)", out).group(1))
    gref = 12683.103 / 1300.0 * (math.pi * 0.5) ** 1.5
    check("Gamma diagnostic matches definition", abs(g / gref - 1) < 1e-4, "%.4f vs %.4f" % (g, gref))


def t_hs_reproduces_thesis(work):
    exe = build_static(work, "monodisperse", {"SYSTEM": 1, "APPROXIMATION": 1}, "mono_hs")
    cwd = os.path.join(work, "v1")
    os.makedirs(cwd)
    run([exe, "0.71", "20", "1024"] + BASIS_HS, cwd, {"LM_PICMAX": "300"})
    mine = read2(latest("Sq__*.dat", cwd))
    ref = read2(os.path.join(DATA, "thesis_hs_mono_py_0.71_Sq.dat"))
    d = max(abs(a[1] - b[1]) for a, b in zip(mine, ref))
    check("V1 HS-PY n=0.71 reproduces thesis S(q)", d <= 1e-8, "max|dS|=%.2e" % d)


def nep(exe, cwd, sqfile, eta_or_p, is_eta=True):
    os.makedirs(cwd)
    name = os.path.basename(sqfile)
    shutil.copy(sqfile, os.path.join(cwd, name))
    p = 4.0 * eta_or_p / math.pi if is_eta else eta_or_p
    n = sum(1 for l in open(sqfile) if l.strip())
    rc, out = run([exe, "500", repr(p), "50", name, str(n)], cwd)
    return rc, read2(os.path.join(cwd, "f" + name[1:]))


def t_nep(work):
    q0 = build_nep(work, 0)
    q1 = build_nep(work, 1)
    rc, f = nep(q0, os.path.join(work, "v4"), os.path.join(DATA, "thesis_sq_0.9"), 0.9, is_eta=False)
    ref = read2(os.path.join(DATA, "thesis_fq_0.9"))
    d = max(abs(a[1] - b[1]) for a, b in zip(f, ref))
    check("V4 QUAD=0 reproduces thesis fq_0.9", d <= 1e-6, "max|df|=%.2e" % d)
    rc, f = nep(q1, os.path.join(work, "v5a"), os.path.join(DATA, "py_sq_eta0.691"), 0.691)
    check("V5 QUAD=1: eta=0.691 is liquid", max(x[1] for x in f) < 1e-6 and rc != 0, "max f=%.2e" % max(x[1] for x in f))
    rc, f = nep(q1, os.path.join(work, "v5b"), os.path.join(DATA, "py_sq_eta0.692"), 0.692)
    check("V5 QUAD=1: eta=0.692 is glass", max(x[1] for x in f) > 0.5 and rc == 0, "max f=%.3f" % max(x[1] for x in f))


def main():
    full = "--full" in sys.argv
    work = tempfile.mkdtemp(prefix="lm_tests_")
    t_binary_compiles(work)
    t_dipolar_calibration_and_no_false_convergence(work)
    if full:
        t_hs_reproduces_thesis(work)
        t_nep(work)
    failed = [n for n, ok in RESULTS if not ok]
    print("\n%d/%d passed%s" % (len(RESULTS) - len(failed), len(RESULTS), "" if full else " (fast; --full for V1/V4/V5)"))
    shutil.rmtree(work, ignore_errors=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
