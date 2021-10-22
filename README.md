
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jmmm

<!-- badges: start -->
<!-- badges: end -->

The goal of jmmm is to fit multivariate mixed models and apply a
discriminant analysis frame work to predict events from the longitudinal
evolution of biomarker data.

The package was adapted from the work by Fieuws *et al.* [Fieuws
Biostatistics
2008](https://doi-org.ru.idm.oclc.org/10.1093/biostatistics/kxm041) In
summary, we fit a high-dimensional joint model by first fitting all
possible pairwise combinations of univariate mixed models within a
subset of patients with and without graft failure. Note that the subset
of subjects without the event of interest should complete the entire
follow-up, and thus no longer be at risk of the event. Subjects who did
not have the event of interest, but who have not yet accrued full
follow-up, should be excluded from model training. The pairwise mixed
models were fitted using the mixed\_model() function from the
[GLMMAdaptive package (version 0.8) developed by
Rizopoulos](https://drizopoulos.github.io/GLMMadaptive/index.html) We
join the pairwise models into a multivariate mixed model by averaging
the parameter estimates. Finally, we join the multivariate mixed models
to the event of interest using a discriminant analysis approach.
Breaking up the large multivariate model into pairwise combinations of
univariate models reduces computational burden. However, as a result of
the many possible pairwise combinations of biomarker evolutions, some
parameters for each longitudinal marker were estimated multiple times.
Therefore, estimates from the pairwise models are averaged to obtain
single estimates for all parameters in the multivariate mixed model. In
addition, subject level contributions to the first and second
derivatives for the pairwise models are averaged and used to generate a
sandwich estimator of the standard errors for the parameters. The
multivatiate mixed models for the event and survivors are used to
construct likelihood profiles for marker evolutions. In turn, the
profiles are used to generate predictions by applying Bayes’ rule.

The work was supported by a grant from the [Dutch Kidney
Foundation](https://nierstichting.nl/) - Grant 19OK003.

## Installation

You can install the development version of `jmmm` from github using the
devtools package:

`devtools::install_github("JanvandenBrand/jmmm")`

## Example

A basic example with simulated data is included below.

``` r
library(jmmm)
# The next libraries are to run the example
library(MASS)
library(gridExtra)
library(tidyverse)
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.3     v dplyr   1.0.7
#> v tidyr   1.1.4     v stringr 1.4.0
#> v readr   2.0.2     v forcats 0.5.1
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::combine() masks gridExtra::combine()
#> x dplyr::filter()  masks stats::filter()
#> x dplyr::lag()     masks stats::lag()
#> x dplyr::select()  masks MASS::select()
```

### Simulate data

First, generate some play data.

``` r
set.seed(1234)
n <- 100 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 5 # maximum follow-up time

# we generate a data frame with the following design: 
# everyone has a baseline measurment, and then measurements at random follow-up times
df <- data.frame(id = as.character(rep(seq_len(n), each = K)),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))
# design matrices for the fixed effects
X1 <- model.matrix(~ sex + time, data = df)
X2 <- model.matrix(~ time, data = df)
X3 <- model.matrix(~ sex + time, data = df)
X4 <- model.matrix(~ time, data = df)
# Random intercepts model
Z1 <- 
  Z2 <- 
  Z3 <- 
  Z4 <- model.matrix(~ 1, data = df)
# fixed effects coefficients 
betas1 <- c(-0.67, -0.25, 0.24) 
betas2 <- c(-1.15, 0.3)  
betas3 <- c(1.34, -0.13, -0.04)
betas4 <- c(0.25, -0.17)
# random effects matrix
D <- c(0.8, 0.7, 0.6, 0.5)
resid <- c(0.5, 0.87) # standard deviation for error terms for y1 and y2 

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D[1])), rnorm(n, sd = sqrt(D[2])), rnorm(n, sd = sqrt(D[3])), rnorm(n, D[4]))
eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b[as.integer(df$id), 1, drop = FALSE]))
eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b[as.integer(df$id), 2, drop = FALSE])) 
eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b[as.integer(df$id), 3, drop = FALSE])) 
eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b[as.integer(df$id), 4, drop = FALSE])) 

# simulate normal longitudinal data
df$y1 <- rnorm(n * K, mean = eta_y1, sd = resid[1])
df$y2 <- rnorm(n * K, mean = eta_y2, sd = resid[2])
# we set binary outcome from the logistic regression
df$y3 <- as.numeric(as.logical(rbinom(n * K, size = 1, prob = plogis(eta_y3))))
df$y4 <- as.numeric(as.logical(rbinom(n * K, size = 1, prob = plogis(eta_y4))))
# Set failures
df$failure <- rnorm(n * K, mean = -0.5*eta_y1+0.8*eta_y2-0.33*eta_y3+0.15*eta_y4, sd = 1.5*mean(resid)) 
df$failure <- exp(df$failure) / (1 + exp(df$failure))
df$failure <- ifelse(df$failure >= 0.6 & df$time >0, 1, 0)
# Set failure times
subjects <- split(df, df$id)
for (i in seq_along(subjects)) {
  if ( dplyr::summarise(subjects[[i]], max(failure)) == 0 ) {
    subjects[[i]]$time_failure <- unlist(
      dplyr::summarise(subjects[[i]], max(time))
      )["max(time)"]
  } else {
    subjects[[i]]$time_failure <- unlist(
      dplyr::summarise(subjects[[i]] %>% 
                         group_by(failure) %>% 
                         dplyr::filter(failure == 1), min(time))
      )["min(time)"]
  }
  subjects[[i]]$failure <- ifelse(max(subjects[[i]]$failure) == 1, 1, 0)
  subjects[[i]] <- subjects[[i]] %>% dplyr::filter(time <= time_failure)
}
df <- do.call(rbind, subjects)
df <- df %>% mutate(sex = as.numeric(sex))
#scale the parameters
df <- df %>% mutate(y1 = scale(y1),
                    y2 = scale(y2))
```

### Initialization steps

To be able to fit the multivariate mixed model,
GLMMapative::mixed\_model() is used. However, this only accepts a single
outcome vector, not a matrix of outcomes. Therefore, we split the
outcome matrix into pairwise combinations of the longitudinal data.

``` r
pairs <- make_pairs(outcomes = c("y1", "y2", "y3", "y4"))
```

We use a support function to check for the data types of the outcomes.
The result will be used to aid in the call to fit the multivariate mixed
model later on.

``` r
model_info <- test_input_datatypes(data = df, pairs = pairs)
```

Next, we unroll the outcomes matrices of the pairwise combinations of
markers into a single outcome vector.

``` r
df_fail <- df %>% filter(failure == 1)
df_fail_stacked <- stack_data(data = df_fail, id = "id", pairs = pairs, covars = c("time", "sex", "time_failure"))
df_nofail <- df %>% filter(failure == 0)
df_nofail_stacked <- stack_data(data = df_nofail, id = "id", pairs = pairs, covars = c("time", "sex"))
```

Another support function can be used to construct the appropriate model
formulas.

``` r
fixed_nofail <- make_fixed_formula(covars = c("time", "sex"))
fixed_fail <- make_fixed_formula(covars = c("time", "sex", "time_failure"))
random <- make_random_formula(id = "id") 
```

### Fitting the Multivatiate Mixed Models

The actual model fitting is largely done by
GLMMadaptive::mixed\_model(). The mmm\_model() function wraps around it.
First we fit the model for the survivors. Note that the fitting
algorithm uses `future.apply` to paralellize the process in a flexible
manner.

``` r
future::plan("multisession", workers = 4)
```

The `progressr` package has been implemented to provide status updates.
You can use the following to show interactive progress updates.

    progressr::handlers(global=TRUE)
    progressr::handlers(list(
      progressr::handler_progress(
        format=":current/:total (:message) [:bar]",
        width=100,
        complete="+"
      )))

``` r
model_nofail <- mmm_model(fixed = fixed_nofail,
                          random = random,
                          id = "id",
                          data = df_nofail,
                          stacked_data = df_nofail_stacked,
                          pairs = pairs,
                          model_families = model_info,
                          iter_EM = 100,
                          iter_qN_outer = 30,
                          nAGQ = 11)
#> Fitting pairwise models:
#> Retrieving derivatives
#> Compiling model output.
```

``` r
model_fail <- mmm_model(fixed = fixed_fail,
                        random = random,
                        id = "id",
                        data = df_fail,
                        stacked_data = df_fail_stacked,
                        pairs = pairs,
                        model_families = model_info,
                        iter_EM = 100,
                        iter_qN_outer = 30,
                        nAGQ = 11)
#> Fitting pairwise models:
#> Retrieving derivatives
#> Compiling model output.
#> The Variance-Covariance matrix was not a positive definite.
#>             Projecting to the closest positive definite.
#>             Note that your results may be unexpected.
```

### Generate predictions

After the multivariate mixed model has been fitted the
mmm\_predictions() function is used to generate predicted probabilities
of the event of interest at various prediction landmarks. First we need
to generate priors though. In addition, the helper function
`get_outcome_type` is used to be able to generate appropriate
predictions and plots.

``` r
outcomes <- c("y1", "y2", "y3", "y4")
mu <- setNames(rep(0, length(outcomes)), outcomes)
re_samples_nofail <- MASS::mvrnorm(1e3,mu = mu, Sigma = model_nofail$vcov)
re_samples_fail <- MASS::mvrnorm(1e3,mu = mu, Sigma = model_fail$vcov)
prior <- get_priors(data = df,
                    time_failure = "time_failure",
                    failure = "failure",
                    horizon = 5,
                    interval = 1/4)
outcome_types <- get_outcome_type(data = df, 
                                  outcomes)
```

Generate the predictions. Note that the formulas for the fixed and
random effects need to include the longitudinal markers as well as
baseline covariates. In order to get predictions at various times, we
simply apply the function across several landmarks. The
`mmm_predictions` function support parallel processing through the
future.apply package.

``` r
predictions <- mmm_predictions(data=df,
                               outcomes=outcome_types,
                               fixed_formula_nofail="~ y1 + y2 + y3 + y4 + time + sex",
                               random_formula_nofail="~ 1| id",
                               random_effects_nofail=re_samples_nofail,
                               parameters_nofail=model_nofail$estimates,
                               fixed_formula_fail="~ y1 + y2 + y3 + y4 + time + sex + time_failure",
                               random_formula_fail="~ 1 | id",
                               random_effects_fail=re_samples_fail,
                               parameters_fail=model_fail$estimates,
                               time="time",
                               failure="failure",
                               failure_time="time_failure",
                               prior=prior,
                               id="id",
                               landmark=1,
                               horizon=5,
                               interval=1/4)
```

### Plot predictions

The following example shows the predictions for a single subject. This
helps communicate forecasts of the probability of an event conditional
on the evolution of longitudinal data. At present the function can only
be used to plot data for a single subject. However it can generate a
list of grobs for predictions generated at different landmark times.

*The code below works in an interactive R session, but not with
markdown. Markdown makes a clean environment when knitting. I haven’t
been able to figure out why this messes with the grid.arrange()*

    plot <- plot_predictions(predictions=predictions,
                     outcomes=outcomes,
                     id="id",
                     subject="1", 
                     time="time"
    )
    grid.arrange(plot)

&lt;Add later: ROC-AUC and calibration plots.&gt;
