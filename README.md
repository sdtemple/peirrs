# Pair-based Estimators of Infection and Removal Rates

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)

<img src="icon.png" align="center" width="600px"/>

Requires `devtools::install_github("sdtemple/pblas")` from R terminal.

### Next steps
- [] Perform some unit tests of code correctness
- [] Use AI to improve docstrings
- [x] Simulate a lot of data
    - [] Redo the 2 beta, 2 gamma simulation (tricky with partial data for pbla)
    - [] Rerun with median / mean approaches
- [x] Simulator with different $\beta$ and $\gamma$
- [x] Qualitatively check correctness
- [x] Multi type implementation of PBLA
- [x] Multi type implementation of tau-based estimator
- [x] Median imputation in tau moments when possible
- [x] Parametric bootstrap
- [x] Multi type MLE with complete data
- [] Spatial tau-based and PBLA estimators
- [] Spatial simulator
- [] Fixed exposure period implementation of tau-based estimator
- [] Truncated exponential for online inference
- [] Multi type Bayesian with complete data
- [] DAMCMC approaches
- [] Python implementation?
- [] EINNs?