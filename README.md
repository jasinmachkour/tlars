---
output:
  html_document:
    variant: markdown_github
    keep_md: yes
  md_document:
    variant: markdown_github
  pdf_document: default
---

<!-- README.md is generated from README.Rmd. Please edit that file -->




# tlars
**Title**: The Terminated-LARS (T_LARS) algorithm for high-dimensional early terminated forward variable selection

**Description**: It computes the solution path of the Terminated-LARS (T-LARS) algorithm. The T_LARS algorithm is a major building block of the Terminating-Knockoff (T-Knock) Filter. The T-Knock filter paper and the corresponding R package are available at [T-Knock paper](https://arxiv.org/abs/2110.06048) and [R package "tknock"](https://github.com/jasinmachkour/tknock), respectively.


## Installation
You can install the tlars package from [GitHub](https://github.com/jasinmachkour/tlars) with 

``` r
install.packages("devtools")
devtools::install_github("jasinmachkour/tlars")
```

You can open the help pages with

```r
library(tlars)
help(package = "tlars")
?tlars
?tlars_model
```

To cite the package ‘tlars’ in publications use:

```r
citation("tlars")
```


## Quick Start
In the following, we illustrate the basic usage of the tlars package to perform variable selection in sparse and high-dimensional regression settings using the T-LARS algorithm.

**First**, we generate a high-dimensional Gaussian data set with sparse support:

```r
library(tlars)

# Setup
n = 100 # number of observations
p = 300 # number of variables
num_act = 10 # number of true active variables
beta = c(rep(0.75, times = num_act), rep(0, times = p - num_act)) # coefficient vector
actives = which(beta > 0) # indices of true active variables
L = p # number of knockoffs

# Generate Gaussian data
set.seed(123)
X = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y = X %*% beta + stats::rnorm(n)
```

**Second**, we generate a knockoff matrix (null/non-active/dummy predictors) with n rows and L knockoff predictors and append it to the original predictor matrix:

```r
set.seed(1234)
knocks = matrix(stats::rnorm(n * L), nrow = n, ncol = L)
XK = cbind(X, knocks)
```

**Third**, we generate an object of the class tlarsCpp and supply the information that the last "L_val" predictors in "XK" are knockoff predictors:

```r
mod_tlars = tlars_model(X = XK, y = y, num_knocks = ncol(knocks))
#> Created an object of class tlars_cpp... 
#> 		 The first 300 predictors in 'XK' are the original predictors and 
#> 		 the last 300 predictors are knockoffs.
mod_tlars
#> C++ object <0x7f9e5506d600> of class 'tlars_cpp' <0x7f9e5709de70>
```

Finally, we perform a tlars-step on "mod_tlars", i.e., the LARS algorithm is run until "T_stop = 1" knockoffs have entered the solution path and stops there.

**An illustrative example**:

******

**T_stop = 1**:


```r
tlars(model = mod_tlars, T_stop = 1, early_stop = TRUE) # performing T-LARS-step on object "larsObj"
#> Executing T-LARS-step by reference...
#> 		 Finished T-LARS-step... 
#> 			 - The results are stored in the C++ object 'mod_tlars'.
#> 			 - New value of T_stop: 1.
#> 			 - Time elaped: 0.002 sec.
beta_hat = mod_tlars$get_beta() # coefficient vector of original variables and knockoffs at the stopping point

eps = .Machine$double.eps # numerical zero

selected_T_1 = which(abs(beta_hat[seq(p)]) > eps) # indices of selected original variables
selected_knock_T_1 = p + which(abs(beta_hat[seq(p + 1, ncol(XK))]) > eps) # indices of selected knockoff variables

FDP_T_1 = length(setdiff(selected_T_1, actives)) / max(1, length(selected_T_1)) # false discovery proportion (FDP)
TPP_T_1 = length(intersect(selected_T_1, actives)) / max(1, length(actives)) # true positive proportion (TPP)

selected_T_1
#> [1]   2   4   6   7   8  10  74 178
selected_knock_T_1
#> [1] 504
FDP_T_1
#> [1] 0.25
TPP_T_1
#> [1] 0.6
```


## Outlook
The [Terminating-Knockoff (T-Knock) filter](https://arxiv.org/abs/2110.06048) from the R package [tknock](https://github.com/jasinmachkour/tknock) is a fast FDR-controlling method for high-dimensional settings that relies on the T-LARS algorithm implemented in the R package [tlars](https://github.com/jasinmachkour/tlars). It allows to control the FDR at a user defined target level. Check it out if you need a fast and FDR-controlling variable/feature selection method for large-scale high-dimensional settings!

## Documentation
For more information and some examples, please check the
[GitHub-vignette](https://github.com/jasinmachkour/tlars/blob/main/vignettes/BasicUsage.Rmd).


## Links
tlars package: [GitHub-tlars](https://github.com/jasinmachkour/tlars).

tknock package: [GitHub-tknock](https://github.com/jasinmachkour/tknock).

T-Knock paper: https://arxiv.org/abs/2110.06048

README file: [GitHub-readme](https://github.com/jasinmachkour/tlars/blob/main/README.md).

Vignette: [GitHub-vignette](https://github.com/jasinmachkour/tlars/blob/main/vignettes/BasicUsage.Rmd).


















