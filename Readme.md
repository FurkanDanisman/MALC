# MALC

**MALC** is our density-estimation algorithm designed for **binned / rounded / interval-grouped** observations.  
This repository provides code to run MALC and compare it against **three kernel-based baseline approaches** used for binned-data density estimation.

---

## Installation

The package is currently available from the GitHub development branch:

```r
install.packages("remotes")
remotes::install_github("FurkanDanisman/DensOLog", ref = "package-format")
```

After installation:

```r
library(MALC)
```

The package is not yet on CRAN, so this will not work until a CRAN release is available:

```r
install.packages("MALC")
```

---

## Basic Usage

The package API uses observed bin frequencies and grid breaks. The main fitting function does not take raw observations.

```r
library(MALC)

counts <- c(2, 6, 12, 20, 12, 6, 2)
grid <- seq(-3.5, 3.5, length.out = length(counts) + 1)

fit <- MALC(counts, grid)
```

Here, `counts[j]` is the observed frequency in the interval
`[grid[j], grid[j + 1])`, with the usual convention that the final bin includes
the right endpoint.

Evaluate the fitted density, CDF, quantiles, and EM-step mean:

```r
x <- seq(min(grid), max(grid), length.out = 200)

dhat <- dmalc(fit, x)
phat <- pmalc(fit, c(-1, 0, 1))
qhat <- qmalc(fit, c(0.25, 0.5, 0.75))
mhat <- mean_malc(fit)
```

Plot the fitted density:

```r
Plot(fit)
```

To fit the smoothed version, set `smooth = TRUE` when fitting. Downstream
functions automatically follow the smoothing choice stored in the fitted object.

```r
fit_s <- MALC(counts, grid, smooth = TRUE)

dmalc(fit_s, x)
pmalc(fit_s, c(-1, 0, 1))
qmalc(fit_s, c(0.25, 0.5, 0.75))
Plot(fit_s)
```

`mean_malc()` returns the EM-step mean estimate and is independent of the
smoothing choice.

---

## Comparison Methods

The package also includes the three comparison methods used in the simulation
study. All use the same `counts, grid` input format.

```r
bk <- BK2002(counts, grid)
bn <- BinnedNP(counts, grid, bandwidth = "boot")
ks <- KernSmooth(counts, grid)

Plot_Kernel(bk)
Plot_Kernel(bn)
Plot_Kernel(ks)
```

Numerical L2 errors can be computed against a known density:

```r
true_pdf <- function(z) dnorm(z)

Eval(fit, true_pdf)
Eval(bk, true_pdf)
Eval(bn, true_pdf)
Eval(ks, true_pdf)
```

---

## What’s in this repo

1. MALC and three baseline (KernSmooth, binnednp, NKDEBD) algorithms.
2. Simulation codes and results to compare methods.
3. Code to compare methods visually.
4. Application of MALC on real data and synthetic data. 
5. A package-style manual in `Manual/MALC-manual.pdf`, with LaTeX source in `Manual/MALC-manual.tex`.

### Repository structure

- `Algorithms/`  
  Algorithms of **MALC** and the **3 comparison methods**, plus shared utilities (e.g., binning helpers, normalizations, evaluation metrics).

- `Simulations/`  
  Scripts to reproduce simulation studies: data generation, repeated runs, saving results, and summary tables.

- `Plot_Comparison/`  
  Plotting code used to generate comparison figures (overlayed estimated densities, error curves, etc.).

- `Application_Example/`  
  End-to-end example on a real (or illustrative) dataset: data prep → estimation → plots.
  
---

## Methods compared

This repo compares **MALC** with the following methods:

### 1) `KernSmooth` (binned kernel density estimation)
We use functions from the **KernSmooth** R package, which implements kernel smoothing and density estimation routines supporting the framework in *Wand & Jones (1995)*.

**Typical usage in this repo:** binned KDE via `bkde()` (and related helpers).

- Package: `KernSmooth` (Ripley; supporting Wand & Jones)
- Reference (book): Wand & Jones (1995)

---

### 2) `binnednp` (kernel estimation for interval-grouped data)
We use routines from the **binnednp** package, which provides kernel density / distribution estimation and bandwidth selection methods for **interval-grouped data**.

- Package: `binnednp`
- Reference: Barreiro-Ures et al. (2019)

> Note: `binnednp` has been archived/removed from CRAN in recent years; if installation is an issue, you can install from an archive snapshot or use the code paths in this repo that do not require it.

---

### 3) Nonlinear kernel density estimation for binned data (Bernoulli)
We include/replicate the method from:

- **Blower & Kelsall (2002)**, *Bernoulli*: “Nonlinear kernel density estimation for binned data: convergence in entropy.”

Bandwidth selection for this method uses selectors available in base R (`stats`), via the `bandwidth` documentation (e.g., `bw.nrd`, `bw.SJ`, etc.), which references classic selectors such as Silverman’s rule-of-thumb and Sheather–Jones.

---

## References

### KernSmooth

	- Wand, M. P., & Jones, M. C. (1995). Kernel Smoothing. Chapman & Hall.
	- Ripley, B. D. (maintainer). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995). R package.

###  binnednp
	- Barreiro-Ures, D., Francisco-Fernández, M., Cao, R., Fraguela, B. B., Doallo, R., González-Andújar, J. L., & Reyes, M. (2019). Analysis of interval-grouped data in weed science: The binnednp Rcpp package. (Preprint / and published version in Ecology & Evolution).

###  Nonlinear KDE for binned data
	- Blower, G., & Kelsall, J. E. (2002). Nonlinear kernel density estimation for binned data: convergence in entropy. Bernoulli, 8(4), 423–449.

###  Bandwidth selection (R base stats)
	- R stats documentation: Bandwidth Selectors for Kernel Density Estimation (e.g., bw.nrd, bw.SJ, etc.).
(Includes references such as Silverman (1986), Scott (1992), and Sheather & Jones (1991).)
