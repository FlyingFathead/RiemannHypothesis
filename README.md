# RiemannHypothesis
**Riemann Hypothesis Numeric Scanner & Exploratory Tools**  

This repository contains experimental Python scripts for **indefinitely** scanning the imaginary axis, searching for zeros of the Riemann zeta function zeta(1/2 + i * t). These scripts:

- Break the t-axis into chunks (e.g. width 40),  
- Use a sign-change method to locate zero-crossings of the real part of a special function Xi(1/2 + i * t),  
- Refine each crossing via bisection, and  
- Log the discovered zeros to both the console and a file (`riemann_zeros.log`).

> **Disclaimer**: This is **not** a proof of the Riemann Hypothesis. No numeric approach alone can prove RH for all infinitely many zeros. This code simply **verifies** zeros on the line up to arbitrarily large ranges, building *empirical* evidence.

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
  - [riemann_hypothesis_scanner.py](riemann_hypothesis_scanner.py) uses naive partial-sum expansions for large imaginary parts of s.  
  - [riemann_hypothesis_rsi.py](riemann_hypothesis_rsi.py) includes a toy “Riemann–Siegel–like” approximation for somewhat better performance at large |t|.

---

## Requirements

- **Python 3.7+**  
- **mpmath** (`pip install mpmath`)  
- Any standard Python environment should suffice.

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

3. **Run a script**  
   - For a basic indefinite scanner:
     ```bash
     python riemann_hypothesis_scanner.py
     ```
   - For a rough “Riemann–Siegel” version:
     ```bash
     python riemann_hypothesis_rsi.py
     ```
4. **Observe output**  
   - The console shows intervals of potential zero-crossings and refined zeros.  
   - A file named `riemann_zeros.log` also logs these findings.

5. **Adjust parameters**  
   - Inside the scripts, tweak variables like `mp.mp.prec` (precision), `CHUNK_SIZE`, `STEP_SIZE`, `REFINE_TOL`, etc. to suit your performance vs. accuracy needs.  
   - By default, it runs indefinitely (`EXPANSIONS = None`). If you want only a few expansions, set `EXPANSIONS` to a finite integer.

---

## Scripts Overview

### 1. `riemann_hypothesis_scanner.py`

- Uses `flexible_zeta(s)` which calls `mp.zeta(s)` for moderate imaginary parts of s, and a simple partial-sum approach for large ones.  
- Logs each discovered zero with a “PASS” if |Xi| < 1e-10.  
- **Not** optimized for extremely large |t| but demonstrates the indefinite scanning concept.

### 2. `riemann_hypothesis_rsi.py`

- Similar indefinite scanning approach.  
- Attempts a naive “Riemann–Siegel–like Z” function for large |t|.  
- Still a toy version: real Riemann–Siegel formula is more sophisticated (involving a proper theta(t) function, Gram points, etc.).  
- May handle large |t| a bit better than the naive script but is not a robust HPC solution.

---

## Performance Notes

- For truly **large |t|**, these scripts may become slow or lose numeric stability, even with `mp.mp.prec` set high. Real HPC computations require advanced expansions or dedicated methods (like implementing the full Riemann–Siegel formula).  
- A smaller `STEP_SIZE` improves accuracy but increases computation time.  
- Increase `mp.mp.prec` if you see suspicious results at large |t|.

---

## Advanced Ideas

- **Implement a full Riemann–Siegel** formula: The toy approach in `riemann_hypothesis_rsi.py` is incomplete. A correct formula would handle the theta(t) phase function and corrections meticulously.  
- **Check prime gap or zero correlation**: Instead of just scanning, you could do statistical analysis on zero spacings, compare to random matrix theory, etc.  
- **Distributed/HPC**: If you want to push well beyond t ~ 1e6 or more, consider parallelizing or using advanced algorithms (Odlyzko’s approach, Gram points, etc.).

---

## License

These scripts are here, and if you crack the million-dollar prize, send me at least 50% (contact: `flyingfathead@protonmail.com`). Thanks.

---

## Final Note

Enjoy scanning zeros **ad infinitum**! Just remember:

> **This is not** a proof of RH. It’s *empirical verification* that zeros up to the scanned range appear on the line.

If you discover any issues or have suggestions for improvements, feel free to open an issue or send a pull request. Happy exploring!
