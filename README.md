[![DOI](https://zenodo.org/badge/122840016.svg)](https://zenodo.org/badge/latestdoi/122840016)



# `R` and `TMB` codes for fitting multievent SIR models to capture-recapture data

We provide here the `R` and `TMB` codes associated with the paper *Eco-epidemiological capture-recapture parameter estimates are robust to variation in infection state uncertainty* by Sarah Benhaiem, Lucile Marescot, Heribert Hofer, Marion L. East, J-D Lebreton, Stephanie Kramer-Schadt, Olivier Gimenez.

First, we present a [benchmark analysis](https://github.com/oliviergimenez/sir_multievent/blob/master/benchmarking_MECMR.md) to compare a native `R` vs. `TMB` ([Kristensen et al. 2016](https://www.jstatsoft.org/article/view/v070i05)) implementation of the optimization procedure.

Second, we provide the `R` and `TMB` codes to [run a Monte Carlo simulation study](https://github.com/oliviergimenez/sir_multievent/blob/master/biasandmse_on%20R_with%20recovery%20probability.R) to assess the bias and mean square error of the multievent SIR capture-recapture model parameters estimators.
