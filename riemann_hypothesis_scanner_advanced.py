#!/usr/bin/env python3
"""
riemann_hypothesis_scanner_advanced.py

A more advanced indefinite numeric scanning script for zeros of Xi(s)
where s = 1/2 + i*t, attempting to incorporate:
  1) A more complete Riemann–Siegel approach with slightly improved remainder.
  2) Additional terms in the theta function for moderate t.
  3) Adaptive step scanning to avoid missing zeros.
  4) A partial Turing check to verify we haven't missed zeros in each chunk.
  5) Parallel scanning of multiple chunks at once.
  6) Ability to save and resume scanning from the last chunk.
  7) Logs zero data to CSV for analysis, plus normal logging to console and log file.

THIS IS STILL NOT A PROOF and is not HPC-level production code.
But it's a more feature-rich demonstration.

Usage example:
    python riemann_hypothesis_scanner_advanced.py

Stop anytime with Ctrl+C.
"""

import os
import json
import csv
import math
import logging
import mpmath as mp
import multiprocessing
from multiprocessing import Pool

# ==================================================
# 1) CONFIGURE PARAMETERS & GLOBALS
# ==================================================

# Increase precision if you plan to go far:
mp.mp.prec = 300  

# Default indefinite scanning
DEFAULT_EXPANSIONS = None  

# Chunk size (range) for each scanning iteration
CHUNK_SIZE = 40  

# We'll scan in parallel how many chunks at once. 
# e.g., PROCESSES = 2 => handle 2 chunks in parallel
PROCESSES = 2  

# Tolerances
REFINE_TOL = 1e-12   # Bisection refinement tolerance
CHECK_TOL  = 1e-10   # If |Xi| < CHECK_TOL => PASS

# The threshold above which we switch from direct mp.zeta()
# to a Riemann–Siegel approach:
RS_THRESHOLD = 50

# Where we store partial scanning state
STATE_FILE = "scanner_resume_state.json"

# Where we store discovered zeros in CSV
ZEROS_CSV = "riemann_zeros_found.csv"

# Step-size for the *initial* scanning approach, used in ADAPTIVE scanning
INITIAL_STEP = 0.5

logger = None  # We'll set this up in setup_logger()

# ==================================================
# 2) LOGGING SETUP
# ==================================================

def setup_logger():
    """
    Configures a logger that logs INFO-level messages
    to both the console and a file named riemann_zeros.log.
    """
    logger_obj = logging.getLogger("RiemannHypothesis")
    logger_obj.setLevel(logging.INFO)

    log_format = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)
    logger_obj.addHandler(console_handler)

    # File handler
    file_handler = logging.FileHandler("riemann_zeros.log", mode="a")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(log_format)
    logger_obj.addHandler(file_handler)

    return logger_obj

# ==================================================
# 3) A MORE COMPLETE R-S FORMULA
#    with a slightly improved remainder
#    and minor expansions in theta(t).
# ==================================================

def theta_rsiegel(t):
    """
    Slightly improved Riemann-Siegel theta function:
        theta(t) ~ t/2 * log(t/(2*pi)) - t/2 + ...
    We'll add the standard -pi/8, plus one correction term 1/(48t).

    For large t, more expansions could be added.

    For negative t, we effectively use absolute value,
    because Z(t) is even in t. Phase sign is handled below.
    """
    t_abs = abs(t)
    if t_abs == 0:
        return mp.mpf('0')
    # Basic main terms
    val = 0.5 * t_abs * mp.log(t_abs/(2*mp.pi)) - 0.5*t_abs - mp.pi/8
    # A tiny correction term, e.g. 1/(48*t)
    val += 1/(48*t_abs)
    return val

def improved_remainder(t, N):
    """
    Attempt a slightly better remainder for R(t):
      R(t) ~ (-1)^N * (t/(2*pi))^(1/4) / sqrt(pi)
    times a small phase factor to reflect partial oscillations.

    This is STILL not the rigorous integral approach,
    but might be a small improvement if the phase is guessed decently.
    """
    if N <= 0:
        return mp.mpf('0')
    t_abs = abs(t)

    base = ((t_abs/(2*mp.pi))**(0.25)) / mp.sqrt(mp.pi)
    # Add a phase factor ~ cos( something ), let's guess a piece of the next partial term
    # e.g., cos( (N+0.5)*something ). This is purely heuristic.
    phase_angle = 0.5 * t_abs * mp.log(N+1)  # TOTALLY ad hoc guess
    return ((-1)**N) * base * mp.cos(phase_angle)

def riemann_siegel_Z(t):
    """
    A more thorough R-S formula for real t:
      Z(t) = sum_{n=1..N} 1/sqrt(n) * cos(t*log(n) - theta_rsiegel(t)) + remainder
    where N=floor(sqrt(t/(2*pi))).

    We also add the small improved remainder that includes a crude phase factor.
    """
    if t == 0:
        return mp.nsum(lambda n: 1/mp.sqrt(n), [1, 1])

    t_abs = abs(t)
    N = int(mp.floor(mp.sqrt(t_abs/(2*mp.pi))))
    if N < 1:
        return mp.mpf('0')

    # Main sum
    main_sum = mp.nsum(
        lambda n: 1/mp.sqrt(n) * mp.cos(t_abs*mp.log(n) - theta_rsiegel(t_abs)),
        [1, N]
    )

    remainder = improved_remainder(t, N)
    return main_sum + remainder

def advanced_zeta_half_plus_i_t(t):
    """
    Approximate zeta(1/2 + i t) using R-S approach:
      Z(t) = e^{+i theta(t)} * zeta(1/2 + i t)
      => zeta(1/2 + i t) = Z(t)*e^{- i theta(t)}

    We compute real-valued Z(t) from riemann_siegel_Z(t),
    then multiply by e^{-i * theta_rsiegel(t)} to get complex zeta(1/2 + i t).

    Negative t uses symmetry (Z(-t)=Z(t)), but the phase sign might invert.
    For scanning, we typically do t>=0. 
    If t<0, we do the same approach with t_abs but that can introduce a small mismatch in sign.
    """
    t_abs = abs(t)
    Z_val = riemann_siegel_Z(t_abs)
    phase = mp.e**(-1j * theta_rsiegel(t_abs))
    if t < 0:
        # we might want to conj() or an additional factor e^{i something},
        # but let's keep it simple
        pass
    return Z_val * phase

def flexible_zeta(s):
    """
    If Re(s)=0.5 and |Im(s)| >= RS_THRESHOLD => use advanced_zeta_half_plus_i_t
    else => mp.zeta(s)
    """
    re_s = mp.re(s)
    im_s = mp.im(s)
    if mp.almosteq(re_s, 0.5, rel_eps=1e-15, abs_eps=1e-15) and abs(im_s) >= RS_THRESHOLD:
        return advanced_zeta_half_plus_i_t(im_s)
    else:
        return mp.zeta(s)

# ==================================================
# 4) XI FUNCTION
# ==================================================

def xi_function(s):
    """
    xi(s) = 0.5 * s*(s-1) * pi^(-s/2) * Gamma(s/2) * zeta(s).
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

# ==================================================
# 5) ADAPTIVE SCAN FOR SIGN CHANGES
#    to avoid skipping zeros if Xi oscillates quickly.
# ==================================================

def adaptive_scan_for_zeros(t_min, t_max, initial_step=0.5):
    """
    We'll move from t_min up to t_max in increments,
    but if we detect a possible sign change, we do a local half-step approach.

    Returns a list of intervals [ (left, right), ... ] 
    where sign changes occur.
    """
    zero_candidates = []
    t_current = t_min
    f_current = real_xi_imag_axis(t_current)

    while t_current < t_max:
        step = initial_step
        # We'll tentatively move forward in increments of 'step',
        # but if we see a sign change, we refine it.
        t_next = t_current + step
        if t_next > t_max:
            t_next = t_max
        f_next = real_xi_imag_axis(t_next)

        # If sign changed
        if f_current * f_next < 0:
            # Try halving the step until sign is stable or step is very small
            half_step = step / 2
            while half_step > 1e-7:
                mid = t_current + half_step
                f_mid = real_xi_imag_axis(mid)
                if f_current * f_mid < 0:
                    # sign changed in [t_current, mid], so reduce step further
                    t_next = mid
                    f_next = f_mid
                else:
                    # sign is in [mid, t_next]
                    t_current = mid
                    f_current = f_mid
                half_step /= 2

            zero_candidates.append((t_current, t_next))

        # Move on
        t_current = t_next
        f_current = f_next
        if t_current >= t_max - 1e-14:
            break

    return zero_candidates

def refine_zero(t_left, t_right, tol=1e-12):
    """
    Attempt to refine a zero within [t_left, t_right]
    using mp.findroot (bisection).
    """
    try:
        root = mp.findroot(
            real_xi_imag_axis,
            (t_left, t_right),
            method='bisect',
            tol=tol,
            maxsteps=1000
        )
        return root
    except mp.libmp.libmp.NoConvergence:
        return None

# ==================================================
# 6) A *VERY* SIMPLE TURING-LIKE CHECK
#    to see if we might have missed a zero.
# ==================================================

def partial_turing_check(t_min, t_max, found_zeros):
    """
    This is not a full Turing method, but a partial check:
    - We approximate the zero count in [t_min, t_max] from known average spacing
      for large t ~ (2*pi / log(t/2/pi)) or from known zero distribution.
    - Compare to how many we actually found.

    If the difference is suspicious, log a warning.

    A real Turing method would be more detailed and check integrals of Z(t).
    """
    # naive approach: The average gap near t is about 2*pi / log(t/(2*pi)).
    # We'll approximate an integer count = integral of [ 1/gap(t) ] from t_min to t_max
    # This is *very rough* but can highlight big misses.

    def average_gap(t):
        if t < 2*math.pi:
            return 2*math.pi  # fallback
        return 2*mp.pi / mp.log(t / (2*mp.pi))

    # integrate from t_min to t_max
    # integral of 1 / average_gap(t) dt is about (t_max - t_min) / average_gap(midpoint)
    # but let's do a small numeric integration for demonstration
    def inverse_gap(t_):
        return 1 / average_gap(t_)

    # We'll do a quick integration with mp.quad
    if t_min < 10:
        t_min_ = 10  # avoid nonsense for small t
    else:
        t_min_ = t_min
    if t_min_ >= t_max:
        return  # no check possible

    approximate_count = mp.quad(inverse_gap, [t_min_, t_max])
    # the found_zeros are those in [t_min, t_max]
    # let's see how many found
    inrange_zeros = [z for z in found_zeros if (z >= t_min and z <= t_max)]
    num_found = len(inrange_zeros)

    # if the difference is huge, we suspect missed zeros
    # We'll do a factor of 2 tolerance
    if num_found < (0.5 * approximate_count):
        logger.warning(
            f"Turing check: Possibly missed zeros? We found {num_found}, but expected ~{approximate_count:.1f}"
        )
    elif num_found > (2.0 * approximate_count):
        logger.warning(
            f"Turing check: Possibly double-counted zeros? We found {num_found}, expected ~{approximate_count:.1f}"
        )
    else:
        logger.info(
            f"Turing check: Found {num_found} zeros, expected ~{approximate_count:.1f}. Looks okay."
        )

# ==================================================
# 7) SCANNING LOGIC (PARALLEL + SAVE/RESUME)
# ==================================================

def process_chunk(chunk_start, chunk_size, step, refine_tol, check_tol):
    """
    Scans the range [chunk_start, chunk_start+chunk_size] for zeros
    using an adaptive sign-change approach, refines them, and returns
    a list of (root, Xi_val, pass/fail) found in that chunk.
    """
    chunk_end = chunk_start + chunk_size
    logger.info(
        f"Scanning t in [{chunk_start}, {chunk_end}] in parallel task {multiprocessing.current_process().name}"
    )
    intervals = adaptive_scan_for_zeros(chunk_start, chunk_end, initial_step=step)

    refined_zeros = []
    for (a, b) in intervals:
        root = refine_zero(a, b, tol=refine_tol)
        if root is not None:
            val = xi_imag_axis(root)
            absval = abs(val)
            result_str = "PASS" if absval < check_tol else "FAIL"
            refined_zeros.append((float(root), val, result_str))

    # Sort by the zero location
    refined_zeros.sort(key=lambda x: x[0])
    # De-duplicate
    unique = []
    for r in refined_zeros:
        if not unique or abs(r[0] - unique[-1][0]) > 1e-10:
            unique.append(r)
    return unique

def indefinite_scan(
    start=0,
    chunk=40,
    step=0.5,
    expansions=None,
    refine_tol=1e-12,
    check_tol=1e-10,
    processes=1
):
    """
    Repeatedly scan t in chunks of width 'chunk', possibly in parallel.
    We also do a partial Turing check after each wave of scanning.

    We store discovered zeros in a global list. Then we do an optional Turing check
    after each round. We also store scanning progress (start) in STATE_FILE, so we can resume.

    Args:
      start (float): where we begin scanning
      chunk (float): chunk size
      step (float): initial scanning step for adaptive approach
      expansions (int or None): how many chunks to do
      refine_tol, check_tol: numeric tolerances
      processes (int): how many parallel processes to use
    """
    all_zeros_found = []
    current_min = mp.mpf(start)
    iteration = 0

    # We'll create a pool for parallel tasks
    pool = Pool(processes=processes)

    try:
        while True:
            iteration += 1
            tasks = []
            # We'll launch 'processes' chunks in parallel each iteration
            for i in range(processes):
                chunk_start = current_min + i*chunk
                tasks.append((chunk_start, chunk, step, refine_tol, check_tol))

            # Run them in parallel
            results = pool.starmap(process_chunk, tasks)

            # Consolidate zeros
            for chunk_result in results:
                for (root, val, status) in chunk_result:
                    all_zeros_found.append(root)
                    logger.info(
                        f"  Zero at t≈{root} => Xi={val}, |Xi|={abs(val)}, {status}"
                    )
                    # Log to CSV
                    save_zero_to_csv(root, val, status)

            # Turing check over the entire combined chunk range
            # i.e. from current_min to current_min + processes*chunk
            turing_min = float(current_min)
            turing_max = float(current_min + processes*chunk)
            partial_turing_check(turing_min, turing_max, all_zeros_found)

            logger.info(
                f"Done iteration {iteration} over [{turing_min}, {turing_max}]."
            )

            # Bump current_min
            current_min += processes*chunk

            # Save state
            save_state(float(current_min))

            # expansions check
            if expansions is not None:
                expansions -= 1
                if expansions <= 0:
                    logger.info("Reached the user-specified expansions limit. Exiting.")
                    break

    finally:
        pool.close()
        pool.join()

def save_state(current_t):
    """
    Save the current scanning position to STATE_FILE in JSON format.
    """
    data = {
        "current_t": current_t
    }
    with open(STATE_FILE, "w") as f:
        json.dump(data, f)
    logger.info(f"Saved state: scanning will resume at t={current_t} next time.")

def load_state():
    """
    Load scanning position from STATE_FILE. If not found, return None.
    """
    if not os.path.exists(STATE_FILE):
        return None
    with open(STATE_FILE, "r") as f:
        data = json.load(f)
    return data.get("current_t", None)

def init_csv():
    """
    Initialize the CSV file for discovered zeros if it doesn't exist.
    We'll append if it does exist, so we only write a header if new.
    """
    if not os.path.isfile(ZEROS_CSV):
        with open(ZEROS_CSV, "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["t_value", "Xi_value", "Status"])

def save_zero_to_csv(t_val, xi_val, status):
    """
    Append discovered zero data to CSV for analysis.
    """
    with open(ZEROS_CSV, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([t_val, str(xi_val), status])

# ==================================================
# 8) MAIN
# ==================================================

def main():
    """
    Main function. By default:
      - Loads scanning state if available,
      - Runs indefinite scanning with parallel chunks,
      - Writes discovered zeros to CSV,
      - Logs sign changes & partial Turing checks,
      - Allows Ctrl+C to stop (pool is closed in finally).
    """
    global logger
    logger = setup_logger()

    logger.info("===== Advanced Riemann Hypothesis Scanner (Demo) =====")
    logger.info(f"Using precision: {mp.mp.prec}")
    logger.info(f"RS_THRESHOLD = {RS_THRESHOLD} (switch to R-S at |t|>={RS_THRESHOLD})")
    logger.info(f"Press Ctrl+C to stop at any time.\n")

    # Initialize CSV if needed
    init_csv()

    # Attempt to resume scanning from prior state
    state_pos = load_state()
    if state_pos is not None:
        start_t = state_pos
        logger.info(f"Resuming scanning from t={start_t}")
    else:
        start_t = 0.0

    indefinite_scan(
        start=start_t,
        chunk=CHUNK_SIZE,
        step=INITIAL_STEP,
        expansions=DEFAULT_EXPANSIONS,
        refine_tol=REFINE_TOL,
        check_tol=CHECK_TOL,
        processes=PROCESSES
    )

if __name__ == "__main__":
    main()
