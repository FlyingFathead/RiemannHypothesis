# RiemannHypothesis
**Riemann Hypothesis Numeric Scanner & Exploratory Tools**  

This repository contains experimental Python scripts for indefinitely scanning the imaginary axis, searching for zeros of the Riemann zeta function zeta(1/2 + i*t). These scripts:

- Break the t-axis into chunks (for example, width 40),  
- Use a sign-change method to locate zero-crossings of the real part of a special function Xi(1/2 + i*t),  
- Refine each crossing via bisection, and  
- Log the discovered zeros to both the console and a file (`riemann_zeros.log`).

> **Disclaimer**: This is **not** a proof of the Riemann Hypothesis. No numeric approach alone can prove RH for all infinitely many zeros. This code simply verifies zeros on the line up to arbitrarily large ranges, building empirical evidence.

---

## Contents

1. [Features](#features)  
2. [Requirements](#requirements)  
3. [How to Use](#how-to-use)  
4. [Scripts Overview](#scripts-overview)  
5. [Performance Notes](#performance-notes)  
6. [Advanced Ideas](#advanced-ideas)  
7. [License](#license)  

---

## Features

- **Indefinite scanning**: The scripts can run forever (or until you stop them with Ctrl+C), chunk after chunk.  
- **Logging**: They log to both console and `riemann_zeros.log`, giving sign-change intervals and refined zero positions.  
- **Flexible approach**:
  - [`riemann_hypothesis_scanner.py`](riemann_hypothesis_scanner.py) uses a simple Riemann–Siegel formula for large |t|, providing improved accuracy at higher imaginary parts of s compared to naive partial sums.  
  - [`riemann_hypothesis_scanner_advanced.py`](riemann_hypothesis_scanner_advanced.py) builds on this further with:
    - More correction terms in the Riemann–Siegel remainder,  
    - Adaptive step scanning to avoid missing quick zero crossings,  
    - A rudimentary Turing-like check,  
    - Parallelization for scanning multiple chunks at once,  
    - State saving (resume scanning),  
    - Logging zeros to both the console/log and a CSV file.  
  - [`riemann_hypothesis_rsi.py`](riemann_hypothesis_rsi.py) is an older, more basic “Riemann–Siegel–like” partial-sum approach.

All demonstrate the concept of indefinite scanning for zeros, but for truly large |t| or HPC-level performance, more sophisticated methods (or further refinements) are needed.

---

## Requirements

- **Python 3.7+**  
- **mpmath** (install via `pip install mpmath`)  
- A standard Python environment (such as virtualenv or conda) is recommended but not required.

---

## How to Use

1. **Clone this repo**  
   ```bash
   git clone https://github.com/FlyingFathead/RiemannHypothesis.git
   cd RiemannHypothesis
   ```

2. **Install mpmath**  
   ```bash
   pip install mpmath
   ```
   or
   ```bash
   pip install -r requirements.txt
   ```
   if you maintain a `requirements.txt`.

3. **Run a script**  
   - For the **advanced** indefinite scanner (with partial Turing checks, adaptive step, parallel scanning, etc.):
     ```bash
     python riemann_hypothesis_scanner_advanced.py
     ```
   - For a simpler indefinite scanner that uses a basic Riemann–Siegel approach:
     ```bash
     python riemann_hypothesis_scanner.py
     ```
   - For an older partial-sum style “Riemann–Siegel–like” version:
     ```bash
     python riemann_hypothesis_rsi.py
     ```

4. **Observe output**  
   - The console shows intervals of potential zero-crossings and refined zeros.  
   - A file named `riemann_zeros.log` also logs these findings.  
   - For the advanced script, discovered zeros also go to a CSV file (`riemann_zeros_found.csv`) for further analysis.

5. **Adjust parameters**  
   - Inside the scripts, tweak variables like `mp.mp.prec` (precision), `CHUNK_SIZE`, `STEP_SIZE`, `REFINE_TOL`, etc. to suit your performance vs. accuracy needs.  
   - By default, the scripts run indefinitely (`EXPANSIONS = None`). If you want only a few expansions, set `EXPANSIONS` to a finite integer in the code.

---

## Scripts Overview

### 1. `riemann_hypothesis_scanner_advanced.py`

- Includes improved Riemann–Siegel corrections, a basic partial Turing-like check, parallel chunk scanning, adaptive step sign-change detection, and optional resume.  
- Logs zeros to console/log file and also to `riemann_zeros_found.csv`.  
- More accurate than `riemann_hypothesis_scanner.py` at larger |t| and less likely to miss zeros due to oscillations.  
- Still not a fully rigorous HPC solution, but a more feature-rich demonstration.

### 2. `riemann_hypothesis_scanner.py`

- Implements `flexible_zeta(s)`, switching to a simple Riemann–Siegel approach for |t| >= 50.  
- Logs each discovered zero with a “PASS” if |Xi| < 1e-10.  
- More accurate than naive partial sums, though does not include advanced features like adaptive stepping or Turing checks.

### 3. `riemann_hypothesis_rsi.py`

- An older indefinite scanning approach with a “Riemann–Siegel–like Z function” that is even more naive.  
- Good for demonstration but replaced by the newer scripts for practical scanning to larger |t|.

---

## Performance Notes

- For truly **large |t|**, these scripts may become slow or lose numeric stability, even with high `mp.mp.prec`. Real HPC methods require advanced expansions (for example, a full Riemann–Siegel integral remainder or Odlyzko–Schönhage algorithms).  
- A smaller `STEP_SIZE` helps avoid missing zeros but increases computation time. Adaptive stepping (in `riemann_hypothesis_scanner_advanced.py`) is one approach to mitigate that.  
- Adjusting `mp.mp.prec` higher slows down each function call, so find a balance between accuracy and performance.

---

## Advanced Ideas

- **Full Riemann–Siegel**: The “improved remainder” in these scripts is still heuristic. A complete formula includes more detailed integrals and expansions, plus possibly Gram points to systematically isolate each zero.  
- **Turing Method**: The partial Turing-like check here is basic. A real Turing approach provides a proof that no zeros are missed in a given range, often by integrating certain functions.  
- **Statistical analysis**: Once zeros are found, you could study their spacing and compare to random matrix predictions for zero distributions.  
- **Distributed/HPC**: For scanning above t ~ 10^6, consider parallel algorithms (Odlyzko–Schönhage, etc.) and large compute clusters.

---

## License

These scripts are here, and if you crack the million-dollar prize, send me at least 50% (contact: `flyingfathead@protonmail.com`). Thanks.

---

## Final Note

Enjoy scanning zeros **ad infinitum**! Remember:  
> This is not a proof of RH. It’s empirical verification that zeros up to the scanned range appear on the line.

If you discover any issues or have suggestions for improvements, feel free to open an issue or send a pull request. Happy exploring!