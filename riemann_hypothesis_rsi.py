#!/usr/bin/env python3
"""
riemann_hypothesis_rsi.py

An indefinite numeric scanning script with a *very rough* Riemann-Siegel "Z" approximation.
Still NOT a proof—just attempts to handle large t a bit better than naive partial sums.
"""

import mpmath as mp
import logging
import time
import math

mp.mp.prec = 200  # adjustable precision

def setup_logger():
    logger = logging.getLogger("RiemannHypothesisRS")
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File
    fh = logging.FileHandler("riemann_zeros.log", mode="a")
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger

logger = setup_logger()

# ---------------------------------------------------------
# 1) Riemann-Siegel-like Z function approximation (rough!)
# ---------------------------------------------------------

def approx_riemann_siegel_Z(t):
    """
    A toy approximation of the Riemann-Siegel "Z" function for large t.
    Z(t) ~ 2 * sum_{n=1..floor(sqrt(t/(2*pi)))} cos(t log n) + ...
    This is quite incomplete but shows the idea.
    If you want a real Riemann-Siegel, you'd implement the full "theta(t)" phase, etc.
    """
    # For demonstration, let's just do a partial sum up to N = floor(sqrt(t / (2*pi))).
    # Then correct for the remainder *very* crudely.
    # A real R-S approach is more intricate (the "Gram points," the full phase function, etc.).
    N = int(math.floor(math.sqrt(abs(t)/(2*math.pi))))
    # sum cos(t log n)
    ssum = mp.nsum(lambda n: mp.cos(t * mp.log(n)), [1, N])
    return ssum  # super naive

def flexible_zeta(s):
    """
    For large Im(s), attempt a Riemann-Siegel "Z" approximation.
    Otherwise, fallback on mp.zeta(s).
    """
    t = abs(mp.im(s))
    if t < 100:
        return mp.zeta(s)
    else:
        # We do "Z(t)" => interpret real(s)=1/2 only, ignoring the rest
        # Then approximate. This is not a direct zeta(s) but a stand-in function
        # that *behaves similarly* for large t on Re(s) = 1/2.
        # Real code would require the full Riemann-Siegel formula. 
        # This is just a toy version:
        # zeta(1/2 + i t) ~ Z(t)?  Actually we need to factor out the principal part...
        # We'll just pretend "Z(t)" stands in for the magnitude of zeta(1/2 + i t).
        # Then return a complex number with that magnitude and 0 phase (toy).
        # Real code does more advanced manipulations. 
        magnitude_approx = approx_riemann_siegel_Z(t)
        # We'll guess a tiny imaginary part to mimic fluctuations
        return magnitude_approx + mp.j*(1e-50*magnitude_approx)

def xi_function(s):
    """
    xi(s) = 0.5 * s*(s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
    Using flexible_zeta(s) for large |Im(s)|.
    """
    prefactor = 0.5 * s * (s-1)
    pi_factor = mp.power(mp.pi, -s/2)
    gamma_factor = mp.gamma(s/2)
    z_val = flexible_zeta(s)
    return prefactor * pi_factor * gamma_factor * z_val

def xi_imag_axis(t):
    s = mp.mpf("0.5") + mp.j*t
    return xi_function(s)

def real_xi_imag_axis(t):
    return mp.re(xi_imag_axis(t))

def scan_for_zeros(t_min, t_max, step=0.5):
    zero_candidates = []
    t_current = t_min
    f_current = real_xi_imag_axis(t_current)

    while t_current < t_max:
        t_next = t_current + step
        f_next = real_xi_imag_axis(t_next)
        if f_current * f_next < 0:
            zero_candidates.append((t_current, t_next))

        t_current = t_next
        f_current = f_next
    return zero_candidates

def refine_zero(t_left, t_right, tol=1e-12):
    try:
        return mp.findroot(real_xi_imag_axis, (t_left, t_right), method='bisect', tol=tol, maxsteps=1000)
    except mp.libmp.libmp.NoConvergence:
        return None

def indefinite_scan(
    start=0,
    chunk=40,
    step=0.5,
    expansions=None,
    refine_tol=1e-12,
    check_tol=1e-10
):
    iteration = 0
    current_min = mp.mpf(start)

    while True:
        iteration += 1
        current_max = current_min + chunk
        logger.info(f"=== Iteration {iteration}: Scanning t in [{current_min}, {current_max}] ===")
        intervals = scan_for_zeros(current_min, current_max, step=step)

        if intervals:
            logger.info("Found sign-change intervals:")
            for (a, b) in intervals:
                logger.info(f"  Potential zero in ~[{a}, {b}]")
        else:
            logger.info("No sign-change intervals found in this chunk.")

        refined_zeros = []
        for (a, b) in intervals:
            root = refine_zero(a, b, tol=refine_tol)
            if root is not None:
                refined_zeros.append(root)

        refined_zeros.sort()
        unique_zeros = []
        for r in refined_zeros:
            if not unique_zeros or abs(r - unique_zeros[-1]) > 1e-10:
                unique_zeros.append(r)

        for r in unique_zeros:
            val = xi_imag_axis(r)
            absval = abs(val)
            result_str = "PASS" if absval < check_tol else "FAIL"
            logger.info(f"  Zero at t≈{r} => Xi={val}, |Xi|={absval}, {result_str}")

        logger.info(f"Done iteration {iteration} over [{current_min}, {current_max}].")

        current_min = current_max
        if expansions is not None:
            expansions -= 1
            if expansions <= 0:
                logger.info("Reached the user-specified expansions limit. Exiting.")
                break

def main():
    # Tweak as you wish:
    T_START = 0
    CHUNK = 40
    STEP = 0.1   # smaller step => finer scanning
    EXP = None   # run indefinitely
    REF_TOL = 1e-12
    CHK_TOL = 1e-10

    logger.info("===== Riemann Hypothesis RSI Playground =====")
    logger.info(f"Precision: {mp.mp.prec}")
    logger.info(f"Starting indefinite scan at t={T_START}, chunk={CHUNK}, step={STEP}")
    logger.info(f"Refine tolerance={REF_TOL}, check tolerance={CHK_TOL}")
    logger.info("Press Ctrl+C to stop at any time.")

    indefinite_scan(
        start=T_START,
        chunk=CHUNK,
        step=STEP,
        expansions=EXP,
        refine_tol=REF_TOL,
        check_tol=CHK_TOL
    )

if __name__ == "__main__":
    main()
