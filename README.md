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
  - [`riemann_hypothesis_scanner.py`](riemann_hypothesis_scanner.py) now uses a simple Riemann–Siegel formula for large |t|, helping with better accuracy for high imaginary parts of s.  
  - [`riemann_hypothesis_rsi.py`](riemann_hypothesis_rsi.py) still has a more basic “Riemann–Siegel–like” partial sum.  

Both approaches demonstrate the concept of indefinite scanning for zeros, but for truly large |t| or HPC-level performance, more refined methods are needed.

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
   - For the updated indefinite scanner with a simple Riemann–Siegel approach:
     ```bash
     python riemann_hypothesis_scanner.py
     ```
   - For the older “Riemann–Siegel–like” partial-sum version:
     ```bash
     python riemann_hypothesis_rsi.py
     ```
4. **Observe output**  
   - The console shows intervals of potential zero-crossings and refined zeros.  
   - A file named `riemann_zeros.log` also logs these findings.

5. **Adjust parameters**  
   - Inside the scripts, tweak variables like `mp.mp.prec` (precision), `CHUNK_SIZE`, `STEP_SIZE`, `REFINE_TOL`, etc. to suit your performance vs. accuracy needs.  
   - By default, it runs indefinitely (`EXPANSIONS = None`). If you want only a few expansions, set `EXPANSIONS` to a finite integer in the code.

---

## Scripts Overview

### 1. `riemann_hypothesis_scanner.py`

- Implements `flexible_zeta(s)`, which calls `mp.zeta(s)` for moderate imaginary parts, but switches to a simple Riemann–Siegel approach for |t| >= 50.  
- This approach leverages a minimal Riemann–Siegel formula, which involves:
  - A theta(t) phase function,  
  - A main sum up to floor(sqrt(t / (2*pi))),  
  - A crude remainder term.  
- Logs each discovered zero with a “PASS” if |Xi| < 1e-10.  
- Much more accurate than naive partial sums at larger |t|, but still not HPC-grade.

### 2. `riemann_hypothesis_rsi.py`

- An earlier indefinite scanning approach with a “Riemann–Siegel–like Z function” that is less fleshed out.  
- Good for demonstration but replaced by the improved approach in `riemann_hypothesis_scanner.py` for practical scanning to larger |t|.

---

## Performance Notes

- For truly **large |t|**, these scripts may become slow or lose numeric stability, even with `mp.mp.prec` set high. More advanced expansions (for example, a full Riemann–Siegel or Odlyzko–Schönhage method) are needed for HPC-level scanning.  
- A smaller `STEP_SIZE` reduces the chance of missing zeros but increases computation time.  
- Increase `mp.mp.prec` if you see suspicious results at large |t|. Keep in mind this slows down each function call.

---

## Advanced Ideas

- **Implement a fully correct Riemann–Siegel formula**: The approach here is still partial. A complete formula includes more detailed remainder approximations and often uses Gram points for systematic zero isolation.  
- **Statistical analysis**: Once zeros are found, you could study their spacing and compare to random matrix theory.  
- **Distributed/HPC**: Pushing above t ~ 10^6 or higher typically requires parallel or specialized algorithms.

---

## License

These scripts are here, and if you crack the million-dollar prize, send me at least 50% (contact: `flyingfathead@protonmail.com`). Thanks.

---

## Final Note

Enjoy scanning zeros **ad infinitum**! Remember:  
> This is not a proof of RH. It’s empirical verification that zeros up to the scanned range appear on the line.

If you discover any issues or have suggestions for improvements, feel free to open an issue or send a pull request. Happy exploring!