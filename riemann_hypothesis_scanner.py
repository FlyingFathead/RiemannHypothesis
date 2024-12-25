#!/usr/bin/env python3
"""
riemann_hypothesis_scanner.py

A 'final' indefinite numeric scanning script with logging:
  - Keeps scanning zeros of zeta(1/2 + i*t) chunk by chunk, 
    logging them both to console and a log file.
  - This is NOT a proof. It's an empirical approach.
  - For very large t, consider a more advanced approach (Riemann–Siegel, HPC, etc.).

Usage example:
    python riemann_hypothesis_scanner.py

You can stop anytime with Ctrl+C.
"""

import mpmath as mp
import time
import logging

# =======================
# 1) CONFIGURE PARAMETERS
# =======================

# Increase precision if you plan to go far:
mp.mp.prec = 200  

# By default, let the scanning run indefinitely:
DEFAULT_EXPANSIONS = None  

# How wide each chunk is:
CHUNK_SIZE = 40  

# Step used to look for sign changes:
STEP_SIZE = 0.1  

# Tolerances:
REFINE_TOL = 1e-12  # Bisection refinement tolerance
CHECK_TOL  = 1e-10  # If |Xi| < CHECK_TOL => PASS


# =========================
# 2) LOGGING SETUP FUNCTION
# =========================

def setup_logger():
    """
    Configures a logger that logs INFO-level messages
    to both the console and a file named riemann_zeros.log.
    """
    logger = logging.getLogger("RiemannHypothesis")
    logger.setLevel(logging.INFO)

    log_format = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)

    # File handler
    file_handler = logging.FileHandler("riemann_zeros.log", mode="a")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    return logger

logger = setup_logger()

# ===============================
# 3) SWITCHABLE ZETA APPROACH
# ===============================

def approximate_riemann_siegel_zeta(s):
    """
    A naive partial-sum approach for large imaginary part of s, 
    only for demonstration. Real Riemann-Siegel is more advanced.
    """
    # If Im(s) is not huge, just use mpmath's zeta
    if abs(mp.im(s)) < 50:
        return mp.zeta(s)
    else:
        N = 5000
        return mp.nsum(lambda n: 1/mp.power(n, s), [1, N])

def flexible_zeta(s):
    """
    Decide whether to use direct mp.zeta(s) or approximate partial-sum 
    for large |Im(s)|.
    """
    if abs(mp.im(s)) < 100:
        return mp.zeta(s)
    else:
        return approximate_riemann_siegel_zeta(s)

# ===============================
# 4) XI FUNCTION
# ===============================

def xi_function(s):
    """
    xi(s) = 0.5 * s*(s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
    Here, zeta(s) is replaced by flexible_zeta(s).
    """
    prefactor = 0.5 * s * (s - 1)
    pi_factor = mp.power(mp.pi, -s/2)
    gamma_factor = mp.gamma(s/2)
    z_val = flexible_zeta(s)
    return prefactor * pi_factor * gamma_factor * z_val

def xi_imag_axis(t):
    """Compute Xi(1/2 + i*t)."""
    s = mp.mpf('0.5') + mp.j*t
    return xi_function(s)

def real_xi_imag_axis(t):
    """Real part of Xi(1/2 + i*t)."""
    return mp.re(xi_imag_axis(t))

# ===============================
# 5) ZERO-FINDING UTILITIES
# ===============================

def scan_for_zeros(t_min, t_max, step=0.5):
    """
    Scan [t_min, t_max] in increments of 'step', 
    looking for sign changes in real_xi_imag_axis(t).
    Returns intervals [left, right] where sign changes occur.
    """
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
    """
    Attempt to refine a zero within [t_left, t_right] 
    using mp.findroot (bisection).
    """
    f = real_xi_imag_axis
    try:
        root = mp.findroot(f, (t_left, t_right), method='bisect', tol=tol, maxsteps=1000)
        return root
    except mp.libmp.libmp.NoConvergence:
        return None

# ===============================
# 6) INDEFINITE SCANNING LOGIC
# ===============================

def indefinite_scan(
    start=0,
    chunk=40,
    step=0.5,
    expansions=None,
    refine_tol=1e-12,
    check_tol=1e-10
):
    """
    Keeps scanning intervals of width 'chunk' for zero candidates:
      [start, start+chunk], [start+chunk, start+2*chunk], ...
    If 'expansions' is None, runs indefinitely. 
    Otherwise, stops after 'expansions' intervals.
    """
    iteration = 0
    current_min = mp.mpf(start)

    while True:
        iteration += 1
        current_max = current_min + chunk
        logger.info(
            f"=== Iteration {iteration}: Scanning t in [{current_min}, {current_max}] ==="
        )
        intervals = scan_for_zeros(current_min, current_max, step=step)

        if intervals:
            logger.info("Found sign-change intervals:")
            for (a, b) in intervals:
                logger.info(f"  Potential zero in ~[{a}, {b}]")
        else:
            logger.info("No sign-change intervals found in this chunk.")

        # Refine those intervals
        refined_zeros = []
        for (a, b) in intervals:
            root = refine_zero(a, b, tol=refine_tol)
            if root is not None:
                refined_zeros.append(root)

        # Deduplicate
        refined_zeros.sort()
        unique_zeros = []
        for r in refined_zeros:
            if not unique_zeros or abs(r - unique_zeros[-1]) > 1e-10:
                unique_zeros.append(r)

        # Evaluate Xi at each root
        for r in unique_zeros:
            val = xi_imag_axis(r)
            absval = abs(val)
            result_str = "PASS" if absval < check_tol else "FAIL"
            logger.info(
                f"  Zero at t≈{r} => Xi={val}, |Xi|={absval}, {result_str}"
            )

        logger.info(
            f"Done iteration {iteration} over [{current_min}, {current_max}]."
        )

        # Move the scanning window up
        current_min = current_max

        # Stop if expansions is set and we've used them up
        if expansions is not None:
            expansions -= 1
            if expansions <= 0:
                logger.info("Reached the user-specified expansions limit. Exiting.")
                break

def main():
    """
    The main function. By default, does indefinite scanning (expansions=None).
    """
    # You can tweak any of these as you wish
    T_START = 0
    CHUNK = CHUNK_SIZE
    STEP  = STEP_SIZE
    EXP = DEFAULT_EXPANSIONS  # None => indefinite
    REF_TOL = REFINE_TOL
    CHK_TOL = CHECK_TOL

    logger.info("===== Riemann Hypothesis Solver (Numeric Playground) =====")
    logger.info(f"Using precision: {mp.mp.prec}")
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

