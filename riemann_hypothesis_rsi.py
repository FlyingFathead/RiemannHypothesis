#!/usr/bin/env python3
"""
riemann_hypothesis_rsi.py

An indefinite numeric scanning script with a more precise Riemann-Siegel 'Z' approximation
for large imaginary parts of s = 1/2 + i t.

STILL NOT A PROOF—just attempts to handle large t more accurately than naive partial sums,
by implementing some of the Riemann–Siegel formula structure.

1) For small |t|, we fall back on mpmath's zeta(s).
2) For large |t|, we approximate zeta(1/2 + i t) via a Riemann–Siegel-like sum and a
   simplified remainder, plus the Riemann–Siegel theta(t) phase.

This script scans the imaginary axis, looking for zeros by sign changes in the real part
of Xi(s). Results are logged to console and to a log file.
"""

import mpmath as mp
import logging
import math

# Adjust precision as needed
mp.mp.prec = 200

def setup_logger():
    """
    Creates a logger that prints to console and also appends to 'riemann_zeros.log'.
    """
    logger = logging.getLogger("RiemannHypothesisRS")
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console Handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File Handler
    fh = logging.FileHandler("riemann_zeros.log", mode="a")
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger

logger = setup_logger()

# ---------------------------------------------------------
# 1) Riemann-Siegel formula components
# ---------------------------------------------------------

def theta(t):
    """
    Riemann-Siegel theta function, approximate form:
      theta(t) ~ (t/2)*log(pi*t) - (t/2) - (pi/8) + [smaller corrections].
    This is already decently accurate for large t.
    """
    # You can add more terms for even better accuracy at large t.
    t_abs = abs(t)
    return 0.5 * t_abs * mp.log(mp.pi * t_abs) - 0.5 * t_abs - mp.pi/8

def riemann_siegel_z(t):
    """
    A more precise Riemann-Siegel approximation for:
      Z(t) = e^{i theta(t)} * zeta(1/2 + i t)
    which is real-valued for real t, ignoring small imaginary errors.

    Here we implement the standard Riemann–Siegel sum:
      Z(t) = sum_{n=1 to N} (1/sqrt(n)) * cos(t*log(n) - theta(t)) + R(t),
    where N = floor(sqrt(t/(2*pi))) for t>0, plus a simple remainder approximation R(t).

    NOTE: For negative t, we use Z(-t) = Z(t), i.e. even function symmetry.
    """
    # The function is even, so we reduce to t >= 0:
    t_abs = abs(t)
    if t_abs == 0:
        # Z(0) ~ sum_{n=1 to 0}... => no sum, but let's define Z(0) carefully
        # This is a corner case; you could do direct zeta(1/2) etc. if you like.
        return mp.nsum(lambda n: 1/mp.sqrt(n) * mp.cos(-theta(0)), [1,1])

    # N = floor( sqrt(t/(2*pi)) ), standard R-S partition
    N = int(mp.floor(mp.sqrt(t_abs/(2*mp.pi))))
    # Main sum
    main_sum = mp.nsum(
        lambda n: 1/mp.sqrt(n) * mp.cos(t_abs*mp.log(n) - theta(t_abs)),
        [1, N]
    ) if N >= 1 else mp.mpf('0')

    # Crude remainder (you can refine this if needed).
    # A better expression involves integrals or expansions, but let's do a simple approach:
    #   R(t) ~ (-1)^N * 1/sqrt(pi) * ( (t / 2pi)^(1/4) )
    # Possibly times some small phase factor, but typically it's real for Z(t).
    if N > 0:
        remainder_factor = (t_abs/(2*mp.pi))**(0.25) / mp.sqrt(mp.pi)
        R = ((-1)**N) * remainder_factor
    else:
        R = mp.mpf('0')

    # Construct the approximate Z(t).
    Z_val = main_sum + R
    return Z_val

# ---------------------------------------------------------
# 2) Flexible zeta(s)
# ---------------------------------------------------------

def flexible_zeta(s):
    """
    For large Im(s), approximate zeta(1/2 + i t) via a Riemann-Siegel approach.
    Otherwise, fallback on mp.zeta(s).

    We'll interpret 'large' as |Im(s)| >= 100 by default. You can adjust the threshold.
    """
    t = abs(mp.im(s))
    # If Re(s) != 0.5, we won't even try to use the R-S formula
    # (this code is specialized for the critical line).
    if (mp.re(s) != 0.5) or (t < 100):
        return mp.zeta(s)
    else:
        # We want Z(t) = e^{i theta(t)} zeta(1/2 + i t).
        # So zeta(1/2 + i t) = e^{-i theta(t)} * Z(t).
        Z_t = riemann_siegel_z(t)
        # e^{-i theta(t)} factor:
        # Real t => e^{-i theta(t)} = cos(-theta(t)) + i sin(-theta(t)) = cos(theta(t)) - i sin(theta(t))
        # We'll just do that explicitly:
        phase = mp.e**(-1j * theta(t))
        return Z_t * phase

# ---------------------------------------------------------
# 3) Xi function for scanning
# ---------------------------------------------------------

def xi_function(s):
    """
    Xi(s) = 0.5 * s * (s-1) * (pi^(-s/2)) * Gamma(s/2) * zeta(s).

    We override zeta(s) with flexible_zeta(s) for large |Im(s)|.
    This yields an approximate Xi(s) that is more accurate at large imaginary s.
    """
    prefactor = 0.5 * s * (s - 1)
    pi_factor = mp.power(mp.pi, -s/2)
    gamma_factor = mp.gamma(s/2)
    z_val = flexible_zeta(s)
    return prefactor * pi_factor * gamma_factor * z_val

def xi_imag_axis(t):
    """
    Convenience function: Xi(0.5 + i t).
    """
    s = mp.mpf("0.5") + mp.j*t
    return xi_function(s)

def real_xi_imag_axis(t):
    """
    Real part of Xi(0.5 + i t), used for sign-change scanning.
    """
    return mp.re(xi_imag_axis(t))

# ---------------------------------------------------------
# 4) Zero-scanning logic
# ---------------------------------------------------------

def scan_for_zeros(t_min, t_max, step=0.5):
    """
    Scan the real axis in increments of 'step' for sign changes in real(Xi(0.5 + i t)).
    Return intervals [t, t+step] that show a sign change.
    """
    zero_candidates = []
    t_current = t_min
    f_current = real_xi_imag_axis(t_current)

    while t_current < t_max:
        t_next = t_current + step
        f_next = real_xi_imag_axis(t_next)
        # Check for sign change
        if f_current * f_next < 0:
            zero_candidates.append((t_current, t_next))

        t_current = t_next
        f_current = f_next

    return zero_candidates

def refine_zero(t_left, t_right, tol=1e-12):
    """
    Refine the root within [t_left, t_right] using mpmath.findroot(bisect).
    Returns the zero location (float) or None if no convergence.
    """
    try:
        return mp.findroot(real_xi_imag_axis, (t_left, t_right),
                           method='bisect', tol=tol, maxsteps=1000)
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
    """
    Repeatedly scan t in [current_min, current_min + chunk], searching for sign changes.
    If found, we refine them to get candidate zeros.

    Args:
        start (float): Starting t-value for scanning
        chunk (float): The chunk size for each iteration
        step (float): The scanning step size
        expansions (int or None): Number of chunks to scan before stopping.
            If None, runs indefinitely until a KeyboardInterrupt.
        refine_tol (float): Tolerance for root refinement (bisect method)
        check_tol (float): Used for a 'PASS' vs 'FAIL' check on |Xi(s)| near zero
    """
    iteration = 0
    current_min = mp.mpf(start)

    while True:
        iteration += 1
        current_max = current_min + chunk
        logger.info(f"=== Iteration {iteration}: Scanning t in [{current_min}, {current_max}] ===")

        # Find sign-change intervals
        intervals = scan_for_zeros(current_min, current_max, step=step)

        if intervals:
            logger.info("Found sign-change intervals:")
            for (a, b) in intervals:
                logger.info(f"  Potential zero in ~[{a}, {b}]")
        else:
            logger.info("No sign-change intervals found in this chunk.")

        # Refine each sign-change interval
        refined_zeros = []
        for (a, b) in intervals:
            root = refine_zero(a, b, tol=refine_tol)
            if root is not None:
                refined_zeros.append(root)

        # Sort and deduplicate
        refined_zeros.sort()
        unique_zeros = []
        for r in refined_zeros:
            if (not unique_zeros) or (abs(r - unique_zeros[-1]) > 1e-10):
                unique_zeros.append(r)

        # Report zero results
        for r in unique_zeros:
            val = xi_imag_axis(r)
            absval = abs(val)
            result_str = "PASS" if absval < check_tol else "FAIL"
            logger.info(f"  Zero at t≈{r} => Xi={val}, |Xi|={absval}, {result_str}")

        logger.info(f"Done iteration {iteration} over [{current_min}, {current_max}].\n")

        current_min = current_max
        if expansions is not None:
            expansions -= 1
            if expansions <= 0:
                logger.info("Reached the user-specified expansions limit. Exiting.")
                break

# ---------------------------------------------------------
# 5) Main entrypoint
# ---------------------------------------------------------

def main():
    """
    Main function to kick off the indefinite scanning with configurable parameters.
    """
    # Tweak as you wish:
    T_START = 0        # Start scanning at t=0
    CHUNK = 40         # Each iteration covers a 40-unit range on the imaginary axis
    STEP = 0.1         # Step size for the scanning (smaller = more thorough)
    EXP = None         # How many expansions/chunks to do (None = indefinite)
    REF_TOL = 1e-12    # Tolerance for refining zero location
    CHK_TOL = 1e-10    # Tolerance for checking Xi near zero

    logger.info("===== Riemann Hypothesis RSI Playground =====")
    logger.info(f"Precision: {mp.mp.prec}")
    logger.info(f"Starting indefinite scan at t={T_START}, chunk={CHUNK}, step={STEP}")
    logger.info(f"Refine tolerance={REF_TOL}, check tolerance={CHK_TOL}")
    logger.info("Press Ctrl+C to stop at any time.\n")

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

## // old version; less precise
##
# #!/usr/bin/env python3
# """
# riemann_hypothesis_rsi.py

# An indefinite numeric scanning script with a *very rough* Riemann-Siegel "Z" approximation.
# Still NOT a proof—just attempts to handle large t a bit better than naive partial sums.
# """

# import mpmath as mp
# import logging
# import time
# import math

# mp.mp.prec = 200  # adjustable precision

# def setup_logger():
#     logger = logging.getLogger("RiemannHypothesisRS")
#     logger.setLevel(logging.INFO)

#     formatter = logging.Formatter(
#         "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
#         datefmt="%Y-%m-%d %H:%M:%S",
#     )

#     # Console
#     ch = logging.StreamHandler()
#     ch.setLevel(logging.INFO)
#     ch.setFormatter(formatter)
#     logger.addHandler(ch)

#     # File
#     fh = logging.FileHandler("riemann_zeros.log", mode="a")
#     fh.setLevel(logging.INFO)
#     fh.setFormatter(formatter)
#     logger.addHandler(fh)

#     return logger

# logger = setup_logger()

# # ---------------------------------------------------------
# # 1) Riemann-Siegel-like Z function approximation (rough!)
# # ---------------------------------------------------------

# def approx_riemann_siegel_Z(t):
#     """
#     A toy approximation of the Riemann-Siegel "Z" function for large t.
#     Z(t) ~ 2 * sum_{n=1..floor(sqrt(t/(2*pi)))} cos(t log n) + ...
#     This is quite incomplete but shows the idea.
#     If you want a real Riemann-Siegel, you'd implement the full "theta(t)" phase, etc.
#     """
#     # For demonstration, let's just do a partial sum up to N = floor(sqrt(t / (2*pi))).
#     # Then correct for the remainder *very* crudely.
#     # A real R-S approach is more intricate (the "Gram points," the full phase function, etc.).
#     N = int(math.floor(math.sqrt(abs(t)/(2*math.pi))))
#     # sum cos(t log n)
#     ssum = mp.nsum(lambda n: mp.cos(t * mp.log(n)), [1, N])
#     return ssum  # super naive

# def flexible_zeta(s):
#     """
#     For large Im(s), attempt a Riemann-Siegel "Z" approximation.
#     Otherwise, fallback on mp.zeta(s).
#     """
#     t = abs(mp.im(s))
#     if t < 100:
#         return mp.zeta(s)
#     else:
#         # We do "Z(t)" => interpret real(s)=1/2 only, ignoring the rest
#         # Then approximate. This is not a direct zeta(s) but a stand-in function
#         # that *behaves similarly* for large t on Re(s) = 1/2.
#         # Real code would require the full Riemann-Siegel formula. 
#         # This is just a toy version:
#         # zeta(1/2 + i t) ~ Z(t)?  Actually we need to factor out the principal part...
#         # We'll just pretend "Z(t)" stands in for the magnitude of zeta(1/2 + i t).
#         # Then return a complex number with that magnitude and 0 phase (toy).
#         # Real code does more advanced manipulations. 
#         magnitude_approx = approx_riemann_siegel_Z(t)
#         # We'll guess a tiny imaginary part to mimic fluctuations
#         return magnitude_approx + mp.j*(1e-50*magnitude_approx)

# def xi_function(s):
#     """
#     xi(s) = 0.5 * s*(s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
#     Using flexible_zeta(s) for large |Im(s)|.
#     """
#     prefactor = 0.5 * s * (s-1)
#     pi_factor = mp.power(mp.pi, -s/2)
#     gamma_factor = mp.gamma(s/2)
#     z_val = flexible_zeta(s)
#     return prefactor * pi_factor * gamma_factor * z_val

# def xi_imag_axis(t):
#     s = mp.mpf("0.5") + mp.j*t
#     return xi_function(s)

# def real_xi_imag_axis(t):
#     return mp.re(xi_imag_axis(t))

# def scan_for_zeros(t_min, t_max, step=0.5):
#     zero_candidates = []
#     t_current = t_min
#     f_current = real_xi_imag_axis(t_current)

#     while t_current < t_max:
#         t_next = t_current + step
#         f_next = real_xi_imag_axis(t_next)
#         if f_current * f_next < 0:
#             zero_candidates.append((t_current, t_next))

#         t_current = t_next
#         f_current = f_next
#     return zero_candidates

# def refine_zero(t_left, t_right, tol=1e-12):
#     try:
#         return mp.findroot(real_xi_imag_axis, (t_left, t_right), method='bisect', tol=tol, maxsteps=1000)
#     except mp.libmp.libmp.NoConvergence:
#         return None

# def indefinite_scan(
#     start=0,
#     chunk=40,
#     step=0.5,
#     expansions=None,
#     refine_tol=1e-12,
#     check_tol=1e-10
# ):
#     iteration = 0
#     current_min = mp.mpf(start)

#     while True:
#         iteration += 1
#         current_max = current_min + chunk
#         logger.info(f"=== Iteration {iteration}: Scanning t in [{current_min}, {current_max}] ===")
#         intervals = scan_for_zeros(current_min, current_max, step=step)

#         if intervals:
#             logger.info("Found sign-change intervals:")
#             for (a, b) in intervals:
#                 logger.info(f"  Potential zero in ~[{a}, {b}]")
#         else:
#             logger.info("No sign-change intervals found in this chunk.")

#         refined_zeros = []
#         for (a, b) in intervals:
#             root = refine_zero(a, b, tol=refine_tol)
#             if root is not None:
#                 refined_zeros.append(root)

#         refined_zeros.sort()
#         unique_zeros = []
#         for r in refined_zeros:
#             if not unique_zeros or abs(r - unique_zeros[-1]) > 1e-10:
#                 unique_zeros.append(r)

#         for r in unique_zeros:
#             val = xi_imag_axis(r)
#             absval = abs(val)
#             result_str = "PASS" if absval < check_tol else "FAIL"
#             logger.info(f"  Zero at t≈{r} => Xi={val}, |Xi|={absval}, {result_str}")

#         logger.info(f"Done iteration {iteration} over [{current_min}, {current_max}].")

#         current_min = current_max
#         if expansions is not None:
#             expansions -= 1
#             if expansions <= 0:
#                 logger.info("Reached the user-specified expansions limit. Exiting.")
#                 break

# def main():
#     # Tweak as you wish:
#     T_START = 0
#     CHUNK = 40
#     STEP = 0.1   # smaller step => finer scanning
#     EXP = None   # run indefinitely
#     REF_TOL = 1e-12
#     CHK_TOL = 1e-10

#     logger.info("===== Riemann Hypothesis RSI Playground =====")
#     logger.info(f"Precision: {mp.mp.prec}")
#     logger.info(f"Starting indefinite scan at t={T_START}, chunk={CHUNK}, step={STEP}")
#     logger.info(f"Refine tolerance={REF_TOL}, check tolerance={CHK_TOL}")
#     logger.info("Press Ctrl+C to stop at any time.")

#     indefinite_scan(
#         start=T_START,
#         chunk=CHUNK,
#         step=STEP,
#         expansions=EXP,
#         refine_tol=REF_TOL,
#         check_tol=CHK_TOL
#     )

# if __name__ == "__main__":
#     main()
