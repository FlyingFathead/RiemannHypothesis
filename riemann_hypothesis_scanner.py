#!/usr/bin/env python3
"""
riemann_hypothesis_scanner.py

An indefinite numeric scanning script with logging to find zeros
of Xi(s) where s = 1/2 + i*t, using a more accurate Riemann–Siegel
style approach for large t.

Key Points:
  - This script logs discovered zeros of zeta(1/2 + i*t) chunk by chunk,
    both to console and a log file.
  - It is NOT a proof, just a numeric experiment.
  - For *very* large t, you might refine the Riemann–Siegel formula further,
    or use HPC methods.

Usage example:
    python riemann_hypothesis_scanner.py

Stop anytime with Ctrl+C.
"""

import mpmath as mp
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
# 3) Riemann–Siegel-based Z(t)
# ===============================

def theta_rsiegel(t):
    """
    Riemann-Siegel theta function (basic version):
        theta(t) ~ 0.5 * t * log(pi*t) - 0.5*t - pi/8
    for t > 0. By symmetry, for t < 0, we can take theta(-t) = -theta(t)
    if needed, but we typically just use |t| in Z(t).
    """
    # We use abs(t) to keep it well-defined. For negative t, the function is even anyway.
    t_abs = abs(t)
    return 0.5 * t_abs * mp.log(mp.pi * t_abs) - 0.5 * t_abs - mp.pi/8

def riemann_siegel_Z(t):
    """
    A simplistic Riemann–Siegel Z-function approximation for real t.
    For large t, it is more accurate than naive partial sums.

    Z(t) = sum_{n=1..N} (1 / sqrt(n)) * cos(t*log(n) - theta_rsiegel(t)) + R(t),
    where N = floor(sqrt(t/(2*pi))).

    Remainder R(t) is approximated crudely:
        R(t) ~ (-1)^N * (t/(2*pi))^(1/4) / sqrt(pi)
    """
    if t == 0:
        # Z(0) can be handled separately, but let's just do a small fallback
        return mp.nsum(lambda n: 1/mp.sqrt(n), [1,1])

    t_abs = abs(t)
    # N = floor(sqrt(t/(2*pi)))
    N = int(mp.floor(mp.sqrt(t_abs/(2*mp.pi))))
    if N < 1:
        return mp.mpf('0')

    # Main sum
    main_sum = mp.nsum(
        lambda n: 1/mp.sqrt(n) * mp.cos(t_abs*mp.log(n) - theta_rsiegel(t_abs)),
        [1, N]
    )

    # Crude remainder
    remainder = 0
    if N > 0:
        remainder = ((-1)**N) * (t_abs/(2*mp.pi))**(0.25) / mp.sqrt(mp.pi)

    return main_sum + remainder

def advanced_zeta_half_plus_i_t(t):
    """
    Approximate zeta(1/2 + i t) using the Riemann-Siegel Z(t):
      Z(t) = e^{+i theta(t)} * zeta(1/2 + i t)
    => zeta(1/2 + i t) = Z(t) * e^{-i theta(t)}

    Returns a complex value that matches zeta(1/2 + i t) in magnitude
    (and hopefully roughly in phase) for real t.
    """
    # For negative t, Z(-t)=Z(t), but the phase might invert, so just use absolute value
    t_abs = abs(t)

    # Compute real-valued Z(t_abs)
    Z_val = riemann_siegel_Z(t_abs)

    # Multiply by e^{-i theta(t_abs)}
    # For real t, theta(-t) = -theta(t), but we skip sign issues by using t_abs inside theta
    phase = mp.e**(-1j * theta_rsiegel(t_abs))

    # If t < 0, we'd expect the function to be conjugated, but for scanning zeros along
    # the line, typically we only do t >= 0. If you do negative t, you might consider
    # an extra factor of e^{± i something}, but let's keep it simple.
    return Z_val * phase


# ===============================
# 4) SWITCHABLE ZETA APPROACH
# ===============================

def flexible_zeta(s):
    """
    Decide whether to use direct mp.zeta(s) or a Riemann-Siegel-based
    approximation for s = 1/2 + i t.

    - If Re(s)==0.5 and |Im(s)| >= 50, we call advanced_zeta_half_plus_i_t(t).
    - Otherwise, we just do mp.zeta(s).
    
    You can tweak the threshold to switch from direct zeta() to R-S approach.
    """
    re_s = mp.re(s)
    im_s = mp.im(s)

    if (abs(im_s) >= 50) and (mp.almosteq(re_s, 0.5, rel_eps=1e-15, abs_eps=1e-15)):
        # For large t on the critical line, use R-S
        return advanced_zeta_half_plus_i_t(im_s)
    else:
        # Fallback
        return mp.zeta(s)

# ===============================
# 5) XI FUNCTION
# ===============================

def xi_function(s):
    """
    xi(s) = 0.5 * s * (s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
    We replace zeta(s) with flexible_zeta(s), which uses
    the Riemann–Siegel approach if appropriate.
    """
    prefactor = 0.5 * s * (s - 1)
    pi_factor = mp.power(mp.pi, -s/2)
    gamma_factor = mp.gamma(s/2)
    z_val = flexible_zeta(s)
    return prefactor * pi_factor * gamma_factor * z_val

def xi_imag_axis(t):
    """
    Compute Xi(1/2 + i*t). We do the scanning along Re(s) = 0.5.
    """
    s = mp.mpf('0.5') + mp.j*t
    return xi_function(s)

def real_xi_imag_axis(t):
    """
    Real part of Xi(1/2 + i*t). 
    Used for sign-change detection in scanning.
    """
    return mp.re(xi_imag_axis(t))

# ===============================
# 6) ZERO-FINDING UTILITIES
# ===============================

def scan_for_zeros(t_min, t_max, step=0.5):
    """
    Scan [t_min, t_max] in increments of 'step', 
    looking for sign changes in real_xi_imag_axis(t).
    Returns intervals [left, right] where sign changes occur.

    NOTE: With a fixed 'step', you might miss zeros if
    Xi(s) oscillates quickly. For more robust detection,
    consider smaller steps or an adaptive approach.
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
    using mp.findroot (bisection). Returns the root or None
    if no convergence.
    """
    f = real_xi_imag_axis
    try:
        root = mp.findroot(
            f,
            (t_left, t_right),
            method='bisect',
            tol=tol,
            maxsteps=1000
        )
        return root
    except mp.libmp.libmp.NoConvergence:
        return None

# ===============================
# 7) INDEFINITE SCANNING LOGIC
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

    For each chunk:
      1) We find sign-change intervals.
      2) Refine them with bisection to locate the zero more accurately.
      3) Log results, including a PASS/FAIL if |Xi| < check_tol.
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

        # Deduplicate any close duplicates
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

# ===============================
# 8) MAIN ENTRY POINT
# ===============================

def main():
    """
    The main function. By default, does indefinite scanning (expansions=None).
    
    If you want to start from a higher t, just adjust T_START.
    If you want fewer or more expansions, set EXP to an integer.
    """
    # You can tweak any of these as you wish:
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

## // old version, poorly optimized
# #!/usr/bin/env python3
# """
# riemann_hypothesis_scanner.py

# A 'final' indefinite numeric scanning script with logging:
#   - Keeps scanning zeros of zeta(1/2 + i*t) chunk by chunk, 
#     logging them both to console and a log file.
#   - This is NOT a proof. It's an empirical approach.
#   - For very large t, consider a more advanced approach (Riemann–Siegel, HPC, etc.).

# Usage example:
#     python riemann_hypothesis_scanner.py

# You can stop anytime with Ctrl+C.
# """

# import mpmath as mp
# import time
# import logging

# # =======================
# # 1) CONFIGURE PARAMETERS
# # =======================

# # Increase precision if you plan to go far:
# mp.mp.prec = 200  

# # By default, let the scanning run indefinitely:
# DEFAULT_EXPANSIONS = None  

# # How wide each chunk is:
# CHUNK_SIZE = 40  

# # Step used to look for sign changes:
# STEP_SIZE = 0.1  

# # Tolerances:
# REFINE_TOL = 1e-12  # Bisection refinement tolerance
# CHECK_TOL  = 1e-10  # If |Xi| < CHECK_TOL => PASS


# # =========================
# # 2) LOGGING SETUP FUNCTION
# # =========================

# def setup_logger():
#     """
#     Configures a logger that logs INFO-level messages
#     to both the console and a file named riemann_zeros.log.
#     """
#     logger = logging.getLogger("RiemannHypothesis")
#     logger.setLevel(logging.INFO)

#     log_format = logging.Formatter(
#         "%(asctime)s [%(levelname)s] %(name)s: %(message)s", 
#         datefmt="%Y-%m-%d %H:%M:%S"
#     )

#     # Console handler
#     console_handler = logging.StreamHandler()
#     console_handler.setLevel(logging.INFO)
#     console_handler.setFormatter(log_format)
#     logger.addHandler(console_handler)

#     # File handler
#     file_handler = logging.FileHandler("riemann_zeros.log", mode="a")
#     file_handler.setLevel(logging.INFO)
#     file_handler.setFormatter(log_format)
#     logger.addHandler(file_handler)

#     return logger

# logger = setup_logger()

# # ===============================
# # 3) SWITCHABLE ZETA APPROACH
# # ===============================

# def approximate_riemann_siegel_zeta(s):
#     """
#     A naive partial-sum approach for large imaginary part of s, 
#     only for demonstration. Real Riemann-Siegel is more advanced.
#     """
#     # If Im(s) is not huge, just use mpmath's zeta
#     if abs(mp.im(s)) < 50:
#         return mp.zeta(s)
#     else:
#         N = 5000
#         return mp.nsum(lambda n: 1/mp.power(n, s), [1, N])

# def flexible_zeta(s):
#     """
#     Decide whether to use direct mp.zeta(s) or approximate partial-sum 
#     for large |Im(s)|.
#     """
#     if abs(mp.im(s)) < 100:
#         return mp.zeta(s)
#     else:
#         return approximate_riemann_siegel_zeta(s)

# # ===============================
# # 4) XI FUNCTION
# # ===============================

# def xi_function(s):
#     """
#     xi(s) = 0.5 * s*(s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
#     Here, zeta(s) is replaced by flexible_zeta(s).
#     """
#     prefactor = 0.5 * s * (s - 1)
#     pi_factor = mp.power(mp.pi, -s/2)
#     gamma_factor = mp.gamma(s/2)
#     z_val = flexible_zeta(s)
#     return prefactor * pi_factor * gamma_factor * z_val

# def xi_imag_axis(t):
#     """Compute Xi(1/2 + i*t)."""
#     s = mp.mpf('0.5') + mp.j*t
#     return xi_function(s)

# def real_xi_imag_axis(t):
#     """Real part of Xi(1/2 + i*t)."""
#     return mp.re(xi_imag_axis(t))

# # ===============================
# # 5) ZERO-FINDING UTILITIES
# # ===============================

# def scan_for_zeros(t_min, t_max, step=0.5):
#     """
#     Scan [t_min, t_max] in increments of 'step', 
#     looking for sign changes in real_xi_imag_axis(t).
#     Returns intervals [left, right] where sign changes occur.
#     """
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
#     """
#     Attempt to refine a zero within [t_left, t_right] 
#     using mp.findroot (bisection).
#     """
#     f = real_xi_imag_axis
#     try:
#         root = mp.findroot(f, (t_left, t_right), method='bisect', tol=tol, maxsteps=1000)
#         return root
#     except mp.libmp.libmp.NoConvergence:
#         return None

# # ===============================
# # 6) INDEFINITE SCANNING LOGIC
# # ===============================

# def indefinite_scan(
#     start=0,
#     chunk=40,
#     step=0.5,
#     expansions=None,
#     refine_tol=1e-12,
#     check_tol=1e-10
# ):
#     """
#     Keeps scanning intervals of width 'chunk' for zero candidates:
#       [start, start+chunk], [start+chunk, start+2*chunk], ...
#     If 'expansions' is None, runs indefinitely. 
#     Otherwise, stops after 'expansions' intervals.
#     """
#     iteration = 0
#     current_min = mp.mpf(start)

#     while True:
#         iteration += 1
#         current_max = current_min + chunk
#         logger.info(
#             f"=== Iteration {iteration}: Scanning t in [{current_min}, {current_max}] ==="
#         )
#         intervals = scan_for_zeros(current_min, current_max, step=step)

#         if intervals:
#             logger.info("Found sign-change intervals:")
#             for (a, b) in intervals:
#                 logger.info(f"  Potential zero in ~[{a}, {b}]")
#         else:
#             logger.info("No sign-change intervals found in this chunk.")

#         # Refine those intervals
#         refined_zeros = []
#         for (a, b) in intervals:
#             root = refine_zero(a, b, tol=refine_tol)
#             if root is not None:
#                 refined_zeros.append(root)

#         # Deduplicate
#         refined_zeros.sort()
#         unique_zeros = []
#         for r in refined_zeros:
#             if not unique_zeros or abs(r - unique_zeros[-1]) > 1e-10:
#                 unique_zeros.append(r)

#         # Evaluate Xi at each root
#         for r in unique_zeros:
#             val = xi_imag_axis(r)
#             absval = abs(val)
#             result_str = "PASS" if absval < check_tol else "FAIL"
#             logger.info(
#                 f"  Zero at t≈{r} => Xi={val}, |Xi|={absval}, {result_str}"
#             )

#         logger.info(
#             f"Done iteration {iteration} over [{current_min}, {current_max}]."
#         )

#         # Move the scanning window up
#         current_min = current_max

#         # Stop if expansions is set and we've used them up
#         if expansions is not None:
#             expansions -= 1
#             if expansions <= 0:
#                 logger.info("Reached the user-specified expansions limit. Exiting.")
#                 break

# def main():
#     """
#     The main function. By default, does indefinite scanning (expansions=None).
#     """
#     # You can tweak any of these as you wish
#     T_START = 0
#     CHUNK = CHUNK_SIZE
#     STEP  = STEP_SIZE
#     EXP = DEFAULT_EXPANSIONS  # None => indefinite
#     REF_TOL = REFINE_TOL
#     CHK_TOL = CHECK_TOL

#     logger.info("===== Riemann Hypothesis Solver (Numeric Playground) =====")
#     logger.info(f"Using precision: {mp.mp.prec}")
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

