# rdmrtmb

**This repository is AI generated and under active development. It is not fit for purpose at this time.**

rdmrtmb is an experimental R package for maximum likelihood estimation of the
racing diffusion model of speeded decision making, using automatic
differentiation via [RTMB](https://github.com/kristian-lange/RTMB).

## Status

- AI generated, with limited human review so far.
- Under active development. The API is **not fixed** and will change without
  warning, including function names, arguments, and output structures.
- **No guarantee that the package works as intended.** It has not been
  thoroughly checked, and results should not be trusted for research or
  production use.

Do not use this package for anything you depend on. If you are interested,
watch the repository or [open an issue](https://github.com/singmann/rdmrtmb/issues)
rather than relying on the current state.

## Overview

Every model parameter is described by its own formula, collected with
`rdm_formula()`, and random effects are written inside these formulas with
lme4-style terms such as `(1 | id)`. Random effects are integrated with the
Laplace approximation.

## Installation

```r
# install.packages("remotes")
remotes::install_github("singmann/rdmrtmb")
```

Expect breakage.

## License

GPL-3
