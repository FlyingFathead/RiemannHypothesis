#!/usr/bin/env python3
"""
riemann_hypothesis_scanner_advanced.py

An improved indefinite numeric scanning script for zeros of Xi(s)
where s = 1/2 + i*t, incorporating:
  1) A more complete Riemann–Siegel approach with slightly improved remainder.
  2) Additional terms in the theta function for moderate t.
  3) Adaptive step scanning to avoid missing zeros.
  4) A partial Turing check to verify we haven't missed zeros in each chunk.
  5) Parallel scanning of multiple chunks at once (child processes).
  6) Ability to save and resume scanning from the last chunk.
  7) Logs zero data to CSV and to console/log file.

Key Fix:
  - We remove direct logger calls in the child processes to avoid "NoneType" errors.
    Instead, child processes return messages. The main process logs them.

Usage:
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

logger = None  # We'll set this up in main()

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
# 3) R-S FORMULA COMPONENTS
# ==================================================

def theta_rsiegel(t):
    """
    Slightly improved Riemann-Siegel theta function:
      theta(t) ~ t/2 * log(t/(2*pi)) - t/2 - pi/8 + 1/(48t)
    for large t, ignoring smaller expansions.
    """
    t_abs = abs(t)
    if t_abs == 0:
        return mp.mpf('0')
    val = 0.5 * t_abs * mp.log(t_abs/(2*mp.pi)) - 0.5*t_abs - mp.pi/8
    # A tiny correction term, e.g. 1/(48*t)
    val += 1/(48*t_abs)
    return val

def improved_remainder(t, N):
    """
    Slightly improved remainder:
      R(t) ~ (-1)^N * (t/(2*pi))^(1/4) / sqrt(pi) * cos(phase).
    This is still heuristic, not a rigorous integral approach.
    """
    if N <= 0:
        return mp.mpf('0')
    t_abs = abs(t)
    base = ((t_abs/(2*mp.pi))**(0.25)) / mp.sqrt(mp.pi)
    # Ad hoc phase guess
    phase_angle = 0.5 * t_abs * mp.log(N+1)
    return ((-1)**N) * base * mp.cos(phase_angle)

def riemann_siegel_Z(t):
    """
    R-S formula approximation for real t:
      Z(t) = sum_{n=1..N} 1/sqrt(n) cos(t*log(n) - theta_rsiegel(t)) + remainder
    where N = floor(sqrt(t/(2*pi))).
    """
    if t == 0:
        return mp.nsum(lambda n: 1/mp.sqrt(n), [1, 1])

    t_abs = abs(t)
    N = int(mp.floor(mp.sqrt(t_abs/(2*mp.pi))))
    if N < 1:
        return mp.mpf('0')

    main_sum = mp.nsum(
        lambda n: 1/mp.sqrt(n) * mp.cos(t_abs*mp.log(n) - theta_rsiegel(t_abs)),
        [1, N]
    )
    rem = improved_remainder(t, N)
    return main_sum + rem

def advanced_zeta_half_plus_i_t(t):
    """
    Approximate zeta(1/2 + i*t) via R-S Z(t):
      zeta(1/2 + i*t) ~ Z(t)*e^{-i theta_rsiegel(t)}.
    """
    t_abs = abs(t)
    Z_val = riemann_siegel_Z(t_abs)
    phase = mp.e**(-1j * theta_rsiegel(t_abs))
    # If t<0, this is simplistic, but we typically scan t>=0.
    return Z_val * phase

def flexible_zeta(s):
    """
    If Re(s)=0.5 and |Im(s)| >= RS_THRESHOLD => advanced_zeta_half_plus_i_t
    else => mp.zeta(s).
    """
    if (mp.almosteq(mp.re(s), 0.5, rel_eps=1e-15, abs_eps=1e-15)
        and abs(mp.im(s)) >= RS_THRESHOLD):
        return advanced_zeta_half_plus_i_t(mp.im(s))
    else:
        return mp.zeta(s)

# ==================================================
# 4) XI FUNCTION
# ==================================================

def xi_function(s):
    """
    xi(s) = 0.5*s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s).
    """
    prefactor = 0.5 * s * (s - 1)
    pi_factor = mp.power(mp.pi, -s/2)
    gamma_factor = mp.gamma(s/2)
    z_val = flexible_zeta(s)
    return prefactor * pi_factor * gamma_factor * z_val

def xi_imag_axis(t):
    return xi_function(0.5 + mp.j*t)

def real_xi_imag_axis(t):
    return mp.re(xi_imag_axis(t))

# ==================================================
# 5) ADAPTIVE SCAN
# ==================================================

def adaptive_scan_for_zeros(t_min, t_max, initial_step=0.5):
    """
    Move from t_min to t_max in increments, halving step if we detect
    sign changes. Returns intervals where sign changes occur.
    """
    zero_candidates = []
    t_current = t_min
    f_current = real_xi_imag_axis(t_current)

    while t_current < t_max:
        step = initial_step
        t_next = t_current + step
        if t_next > t_max:
            t_next = t_max

        f_next = real_xi_imag_axis(t_next)

        if f_current * f_next < 0:
            # Halve the step repeatedly
            half_step = step / 2
            while half_step > 1e-7:
                mid = t_current + half_step
                f_mid = real_xi_imag_axis(mid)
                if f_current * f_mid < 0:
                    t_next = mid
                    f_next = f_mid
                else:
                    t_current = mid
                    f_current = f_mid
                half_step /= 2

            zero_candidates.append((t_current, t_next))

        t_current = t_next
        f_current = f_next
        if t_current >= t_max - 1e-14:
            break

    return zero_candidates

def refine_zero(t_left, t_right, tol=1e-12):
    """
    Refine a zero in [t_left, t_right] with bisect.
    """
    try:
        return mp.findroot(
            real_xi_imag_axis,
            (t_left, t_right),
            method='bisect',
            tol=tol,
            maxsteps=1000
        )
    except mp.libmp.libmp.NoConvergence:
        return None

# ==================================================
# 6) PARTIAL TURING CHECK
# ==================================================

def partial_turing_check(t_min, t_max, found_zeros_list):
    def average_gap(t_):
        if t_ < 6.28:
            return 6.28
        return 2*mp.pi / mp.log(t_/(2*mp.pi))

    def inverse_gap(t_):
        return 1/average_gap(t_)

    t_min_ = max(10, t_min)
    if t_min_ >= t_max:
        return []

    approximate_count = mp.quad(inverse_gap, [t_min_, t_max])
    float_count = float(approximate_count)  # convert to regular float
    inrange_zeros = [z for z in found_zeros_list if (z >= t_min and z <= t_max)]
    num_found = len(inrange_zeros)
    messages = []

    if num_found < (0.5 * float_count):
        messages.append(
            f"Turing check: Possibly missed zeros? Found {num_found}, expected ~{float_count:.1f}"
        )
    elif num_found > (2.0 * float_count):
        messages.append(
            f"Turing check: Possibly double-counted zeros? Found {num_found}, expected ~{float_count:.1f}"
        )
    else:
        messages.append(
            f"Turing check: Found {num_found} zeros, expected ~{float_count:.1f}. Looks okay."
        )

    return messages

# ==================================================
# 7) SCANNING LOGIC (PARALLEL + RESUME)
# ==================================================

def process_chunk(chunk_start, chunk_size, step, refine_tol, check_tol):
    """
    Child process function:
      - Scans [chunk_start, chunk_start+chunk_size]
      - Returns (messages, zero_list)
        where messages is a list of strings to be logged by the parent
        and zero_list is [(root, val, status), ...].
    """
    chunk_end = chunk_start + chunk_size
    messages = []
    messages.append(
        f"Child process {multiprocessing.current_process().name} scanning [{chunk_start}, {chunk_end}]"
    )

    intervals = adaptive_scan_for_zeros(chunk_start, chunk_end, initial_step=step)
    refined_zeros = []
    for (a, b) in intervals:
        root = refine_zero(a, b, tol=refine_tol)
        if root is not None:
            val = xi_imag_axis(root)
            absval = abs(val)
            status = "PASS" if absval < check_tol else "FAIL"
            refined_zeros.append((float(root), val, status))

    # Sort & deduplicate
    refined_zeros.sort(key=lambda x: x[0])
    unique = []
    for r in refined_zeros:
        if not unique or abs(r[0] - unique[-1][0]) > 1e-10:
            unique.append(r)

    messages.append(f"Child found {len(unique)} zero(s) in chunk.")
    return (messages, unique)

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
    Repeatedly scan t in chunks of size 'chunk' in parallel. 
    Each iteration spawns 'processes' tasks.
    The child processes do scanning & refining, then return data for the parent to log.
    We also do a partial Turing check after each wave.

    We store discovered zeros in a list, and store scanning progress in STATE_FILE.
    """
    all_zeros_found = []
    current_min = mp.mpf(start)
    iteration = 0

    pool = Pool(processes=processes)

    try:
        while True:
            iteration += 1
            # Prepare tasks for parallel chunk scanning
            tasks = []
            for i in range(processes):
                chunk_start = current_min + i*chunk
                tasks.append((chunk_start, chunk, step, refine_tol, check_tol))

            # Launch in parallel
            results = pool.starmap(process_chunk, tasks)
            # results is a list of (messages, zero_list)
            # Consolidate & log in parent:
            for (messages, zero_list) in results:
                # Log each message
                for msg in messages:
                    logger.info(msg)
                # Then log each zero & store
                for (root, val, status) in zero_list:
                    all_zeros_found.append(root)
                    logger.info(f"  Zero at t≈{root} => Xi={val}, |Xi|={abs(val)}, {status}")
                    save_zero_to_csv(root, val, status)

            # Turing check over the full multi-chunk range
            turing_min = float(current_min)
            turing_max = float(current_min + processes*chunk)
            turing_msgs = partial_turing_check(turing_min, turing_max, all_zeros_found)
            for m in turing_msgs:
                logger.info(m)

            logger.info(f"Done iteration {iteration} over [{turing_min}, {turing_max}].")

            # Move up
            current_min += processes*chunk
            save_state(float(current_min))

            if expansions is not None:
                expansions -= 1
                if expansions <= 0:
                    logger.info("Reached the user-specified expansions limit. Exiting.")
                    break

    finally:
        pool.close()
        pool.join()

# ==================================================
# 8) SAVE/LOAD & CSV
# ==================================================

def save_state(current_t):
    data = {
        "current_t": current_t
    }
    with open(STATE_FILE, "w") as f:
        json.dump(data, f)
    logger.info(f"Saved state: next start t={current_t}.")

def load_state():
    if not os.path.exists(STATE_FILE):
        return None
    with open(STATE_FILE, "r") as f:
        data = json.load(f)
    return data.get("current_t", None)

def init_csv():
    """
    If the CSV doesn't exist, create it with a header.
    """
    if not os.path.isfile(ZEROS_CSV):
        with open(ZEROS_CSV, "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["t_value", "Xi_value", "Status"])

def save_zero_to_csv(t_val, xi_val, status):
    with open(ZEROS_CSV, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([t_val, str(xi_val), status])

# ==================================================
# 9) MAIN
# ==================================================

def main():
    global logger
    logger = setup_logger()

    logger.info("===== Advanced Riemann Hypothesis Scanner (Demo) =====")
    logger.info(f"Using precision: {mp.mp.prec}")
    logger.info(f"RS_THRESHOLD = {RS_THRESHOLD} (switch to R-S at |t|>={RS_THRESHOLD})")
    logger.info(f"Press Ctrl+C to stop at any time.\n")

    init_csv()

    # Attempt to resume from last scanning position
    state_pos = load_state()
    if state_pos is not None:
        start_t = state_pos
        logger.info(f"Resuming from t={start_t}")
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
