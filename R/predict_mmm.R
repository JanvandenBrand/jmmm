# Get predictions from the mmm ---------------------------------------------------

#'  Return a data.frame with the outcome label and type
#'
#' @param data original data
#' @param outcomes a character vector with the names of the longitudinal outcomes
#'
#' @export

get_outcome_type <- function(data, outcomes) {
  out <- lapply(outcomes, function(i) {
    suppressWarnings(
      if (all(sort(unlist(unique(data[paste(i)]))) == c(0, 1), na.rm = TRUE)) {
        outcome_type <- "binary"
      } else if (class(unlist(data[paste(i)])) == "numeric") {
        outcome_type <- "continuous"
      } else {
        stop(paste("The outcome", i, "does not appear to be binary or continuous"))
      }
    )
    data.frame(outcome = i,
               outcome_type = outcome_type)
  })
  do.call(rbind, out)
}

# FIX: grepl for variable names results in poor abstraction and generalization. Parts of variable names may be similar or hold symbols used in regex.
# FIX: models should be different for different longitudinal outcomes

#'  Get Likelihood No Fail
#'
#' Get the parameter estimates from a MMM and combine with data to return contributions to the likelihood at observed time points for non failures
#'
#' @param outcomes a character vector with the labels for the longitudinal outcomes.
#' @param data data.frame with original data in long format
#' @param fixed_formula_nofail a character string with lme4 style formula for the fixed effects
#' @param random_formula_nofail a character string with lme4 style formula for the random effects
#' @param time a character string with the label for the follow-up time variable
#' @param parameters_nofail a data.frame with the parameter labels and estimates returned by mmm_model()
#' @param random_effects_nofail The random effects sampled from a multivariate normal distribution (see MASS::mvrnorm)
#' @param id a character string with subject identifier
#' @param horizon the prediction horizon
#'
#' @import future.apply
#' @import tidyverse
#' @import Matrix
#'
#' @return A list with the likelihood contributions for each subject (rows) at every time point (cols).

get_likelihood_nofail <- function (data,
                                   outcomes,
                                   fixed_formula_nofail = NULL,
                                   random_formula_nofail = NULL,
                                   time,
                                   parameters_nofail,
                                   random_effects_nofail,
                                   id,
                                   horizon,
                                   ...) {
  # cut the fixed and random formulas into parts
  ff <- str_remove(fixed_formula_nofail, "~ ")
  ff <- str_split(ff, " \\+ ") %>% unlist()
  rf <- str_remove(random_formula_nofail, "~ ")
  rf <- str_remove(rf, paste("\\|", id))
  rf <- str_split(rf, "\\+ ") %>% unlist()
  # Select the parameter estimates for fixed effects
  params <- parameters_nofail %>%
    dplyr::filter(grepl("^b", parameter_name)) %>%
    dplyr::select(parameter_name, parameter_estimate)
  # Initialize lists to store the log-likelihoods by outcome for each individual i, at every time point
  # Loop over the patterns. Note that non-failures only have a single pattern as they all survive past the max horizon.
  .s <- horizon
  # Calculate the log-likelihoods for every outcome j
  likelihood <- future_lapply( 1:nrow(unique(data[id])), function(i, ...) {
    # i <- 1
    # Calculate the log-likelihoods for every outcome j
    ll <- lapply(1:nrow(outcomes), function (j) {
      # j <- 1
      # Get the fixed effect parameters
      betas <- params %>%
        dplyr::filter(grepl(paste0("_", outcomes[j, 1], "$"), parameter_name)) %>%
        dplyr::select(parameter_estimate, parameter_name)
      # adjust the order of the parameter estimates to match the fixed_formula argument
      betas <- betas %>%
        mutate(parameter_name = str_replace(parameter_name, "b0_", "intercept"))
      betas <- betas %>%
        mutate(parameter_name = str_remove(parameter_name, paste0("_", outcomes[j, 1], "$")))
      betas <- betas %>%
        mutate(parameter_name = str_remove(parameter_name, paste0("^(..|...)_")))
      betas <- betas[match(c("intercept", ff[-j]), betas$parameter_name), ] %>%
        dplyr::select(parameter_estimate)
      betas <- as.matrix(betas)
      if(!assertthat::noNA(betas)) {
        stop("Not all parameters have matching estimates. Have you correctly specified the the fixed_formula_nofail argument?")
      }
      # Get the residual
      resid <- parameters_nofail %>%
        dplyr::filter(grepl(paste0("resid", paste0("_", outcomes[j, 1],"$")), parameter_name)) %>%
        dplyr::select(parameter_estimate)
      resid <- as.numeric(resid)
      # get the log likelihood for the i-th subject
      # log_likelihood <- vector("list" , length = nrow(unique(data[id])))
      .id <- unlist(unique(data[id]))[i]
      # Get the outcome vector for subject i
      y <- data %>%
        dplyr::filter(all_of(id) == .id & time <= .s) %>%
        dplyr::select(all_of(outcomes[j, 1]))
      y <- as.matrix(y)
      # Using bX assumes similar models for all longitudinal outcomes.
      .ff <- str_remove(ff, paste0("^", outcomes[j, 1], "$"))
      .ff <- .ff[.ff != ""]
      # Get the X matrix for subject i including all data prior to the horizon.
      X <- data %>%
        dplyr::filter(all_of(id) == .id & time <= .s)
      X <- model.matrix(as.formula(paste("~ ", paste(.ff, collapse = " + "))), X)
      # Get the linear predictor for the fixed effect
      xb <- X %*% betas
      # Get the Z matrix for subject i
      .rf <- str_remove(rf, paste0("^",outcomes[j, 1], "$"))
      .rf <- .rf[.rf != ""]
      # Select all data prior to the prediction horizon.
      Z <- data %>%
        dplyr::filter(all_of(id) == .id & time <= .s)
      Z <- model.matrix(as.formula(paste("~", paste(.rf, collapse = " + "))), Z)
      # log_lik is the log likelihood of the j-th outcome for the i-th subject at observed time points
      # for the k-th draw from the random effects.
      if (outcomes$outcome_type[j] == "continuous") {
        # for (k in 1:nrow(random_effect)) {
        log_lik <- lapply(1:nrow(random_effects_nofail), function(k) {
          # get the linear predictor = X*beta + Z*bu
          zb <- rowSums(Z * random_effects_nofail[k, j])
          lp <- xb + zb
          # calculate the log-likelihood
          -log(resid) -0.5 * log(2 * pi) - ( (y - lp)^2 / (2 * resid^2) )
        })
      } else if (outcomes$outcome_type[j] == "binary") {
        log_lik <- lapply(1:nrow(random_effects_nofail), function(k) {
          # lp = X*beta + Z*bu
          zb <- rowSums(Z * random_effects_nofail[k, j])
          lp <- xb + zb
          p <- exp(lp) / (1 + exp(lp))
          lik <- ifelse(y == 0, 1 - p, p)
          lik <- ifelse(lik > 1e-8, lik, 1e-10)
          log(lik)
        })
      }
      # log likelihood for j-th outcome at each observed time point (columns) and all RE draws (rows) for i-th subject
      list("ll" = t(do.call(cbind, log_lik)),
           ".id" = .id,
           "X" = X)
    })
    # return the sum over the outcomes at each observed time point (columns) and all RE  draws (rows) for i-the subject
    sum_ll <- ll[[1]]$ll
    for (j in 2:length(ll)) {
      sum_ll <- sum_ll + ll[[j]]$ll
    }
    # take the exp to obtain likelihood per draw (rows) and per observed time t (columns) for the i-th subject
    l <- exp(sum_ll)
    # take the mean of the likelihoods of each subject over the draws per observed time t
    l <- matrix(colMeans(l)) # this transposes the times to the rows
    data.frame(
      id = ll[[1]]$.id,
      time = ll[[1]]$X[, time],
      likelihoods = l,
      row.names = NULL
    )
  })
  likelihood
}

# FIX: de likelihoods veranderen niet afhankelijk van de failure time...

#' Get Likelihood Fail
#'
#' Get the parameter estimates from a MMM and combine with data to return contributions to the likelihood at observed time points for non failures
#'
#' @param outcomes a character vector with the labels for the longitudinal outcomes.
#' @param data data.frame with original data in long format
#' @param fixed_formula_fail a character string with lme4 style formula for the fixed effects
#' @param random_formula_fail a character string with lme4 style formula for the random effects
#' @param time a character string with the label for the follow-up time variable
#' @param failure_time a character string with the label for the failure time variable
#' @param parameters_fail a data.frame with the parameter labels and estimates returned by mmm_model()
#' @param random_effects_fail The random effects sampled from a multivariate normal distribution (see MASS::mvrnorm)
#' @param id a character string with subject identifier
#' @param horizon the prediction horizon
#' @param interval the intervals at which predictions are to be made
#' @param landmark the landmark time until which data is available.
#'
#' @import future.apply
#' @import tidyverse
#' @import Matrix
#'
#' @return A list with the likelihood contributions for each subject (rows) at every time point (cols).
get_likelihood_fail <- function (data,
                                 outcomes,
                                 fixed_formula_fail = NULL,
                                 random_formula_fail = NULL,
                                 time,
                                 failure_time,
                                 parameters_fail,
                                 random_effects_fail,
                                 id,
                                 landmark,
                                 horizon,
                                 interval,
                                 ...) {
  # Check if there are observations before the landmark for all subjects
  lm <- seq(from = landmark, to = horizon, by = interval)
  # cut the fixed and random formulas into parts
  ff <- str_remove(fixed_formula_fail, "~ ")
  ff <- str_split(ff, " \\+ ") %>% unlist()
  rf <- str_remove(random_formula_fail, "~ ")
  rf <- str_remove(rf, paste("\\|", id))
  rf <- str_split(rf, "\\+ ") %>% unlist()
  # Select the parameter estimates for fixed effects
  params <- parameters_fail %>%
    dplyr::filter(grepl("^b", parameter_name)) %>%
    dplyr::select(parameter_name, parameter_estimate)
  likelihood_pattern <- lapply(lm, function(s, ...) {
    # s <- 1
    # use midpoint rule.
    .s <- (s + 0.5)
    likelihood <- future_lapply(1:nrow(unique(data[id])), function(i, ...) {
      # i <- 1
      # Calculate the log-likelihoods for every outcome j
      ll <- lapply(1:nrow(outcomes), function (j) {
        # j <- 5
        # Get the fixed effect parameters
        betas <- params %>%
          dplyr::filter(grepl(paste0("_", outcomes[j, 1], "$"), parameter_name)) %>%
          dplyr::select(parameter_estimate, parameter_name)
        # adjust the order of the parameter estimates to match the fixed_formula argument
        betas <- betas %>%
          mutate(parameter_name = str_replace(parameter_name, "b0_", "intercept"))
        betas <- betas %>%
          mutate(parameter_name = str_remove(parameter_name, paste0("_", outcomes[j, 1], "$")))
        betas <- betas %>%
          mutate(parameter_name = str_remove(parameter_name, paste0("^(..|...)_")))
        betas <- betas[match(c("intercept", ff[-j]), betas$parameter_name), ] %>%
          dplyr::select(parameter_estimate)
        betas <- as.matrix(betas)
        if(!assertthat::noNA(betas)) {
          stop("Not all parameters have matching estimates. Have you correctly specified the the fixed_formula_nofail argument?")
        }
        # Get the residual
        resid <- parameters_fail %>%
          dplyr::filter(grepl(paste0("resid", paste0("_", outcomes[j, 1],"$")), parameter_name)) %>%
          dplyr::select(parameter_estimate)
        resid <- as.numeric(resid)
        # get the log likelihood for the i-th subject
        # log_likelihood <- vector("list" , length = nrow(unique(data[id])))
        .id <- unlist(unique(data[id]))[i]
        # Get the outcome vector for subject i
        y <- data %>%
          dplyr::filter(all_of(id) == .id & time <= .s) %>%
          dplyr::select(all_of(outcomes[j, 1]))
        y <- as.matrix(y)
        # bX assumes similar models for all longitudinal outcomes.
        .ff <- str_remove(ff, paste0("^", outcomes[j, 1], "$"))
        .ff <- .ff[.ff != ""]
        # Get the X matrix for subject i, select only data prior to the lm time
        X <- data %>%
          dplyr::filter(all_of(id) == .id & time <= .s)
        X <- model.matrix(as.formula(paste("~ ", paste(.ff, collapse = " + "))), X)
        # replace actual failure time with failure time for the prediction
        if (failure_time %in% colnames(X)) {
          X[ , failure_time] <- .s
        }
        # Get the linear predictor for the fixed effect
        xb <- X %*% betas
        # Get the Z matrix for subject i
        .rf <- str_remove(rf, paste0("^",outcomes[j, 1], "$"))
        .rf <- .rf[.rf != ""]
        # Select only data prior to lm time
        Z <- data %>%
          dplyr::filter(all_of(id) == .id & time <= .s)
        Z <- model.matrix(as.formula(paste("~", paste(.rf, collapse = " + "))), Z)
        # replace actual failure time with failure time for the prediction
        if (failure_time %in% colnames(Z)) {
          Z[ , failure_time] <- .s
        }
        # log_lik is the log likelihood of the j-th outcome for the i-th subject at observed time points
        # for the k-th draw from the random effects.
        if (outcomes$outcome_type[j] == "continuous") {
          # for (k in 1:nrow(random_effect)) {
          log_lik <- lapply(1:nrow(random_effects_fail), function(k) {
            # get the linear predictor = X*beta + Z*bu
            zb <- rowSums(Z * random_effects_fail[k, j])
            lp <- xb + zb
            # calculate the log-likelihood
            -log(resid) -0.5 * log(2 * pi) - ( (y - lp)^2 / (2 * resid^2) )
          })
        } else if (outcomes$outcome_type[j] == "binary") {
          log_lik <- lapply(1:nrow(random_effects_fail), function(k) {
            # lp = X*beta + Z*bu
            zb <- rowSums(Z * random_effects_fail[k, j])
            lp <- xb + zb
            p <- exp(lp) / (1 + exp(lp))
            lik <- ifelse(y == 0, 1 - p, p)
            lik <- ifelse(lik > 1e-8, lik, 1e-10)
            log(lik)
          })
        }
        # log likelihood for j-th outcome at each observed time point (columns) and all RE draws (rows) for i-th subject
        list("ll" = t(do.call(cbind, log_lik)),
             ".id" = .id,
             "X" = X)
      })
      # return the sum over the outcomes at each observed time point (columns) and all RE  draws (rows) for i-the subject
      sum_ll <- ll[[1]]$ll
      for (j in 2:length(ll)) {
        sum_ll <- sum_ll + ll[[j]]$ll
      }
      # take the exp to obtain likelihood per draw (rows) and per observed time t (columns) for the i-th subject
      l <- exp(sum_ll)
      # take the mean of the likelihoods of each over the draws per observed time t
      l <- matrix(colMeans(l))
      data.frame(
        id = ll[[1]]$.id,
        time = ll[[1]]$X[, time],
        likelihoods = l,
        pattern = s,
        row.names = NULL
      )
    })
    likelihood
  })
  likelihood_pattern
}

#' Get likelihood profiles for non failures
#'
#' Generate a step function of the likelihood at interval times.
#'
#' @param likelihood_nofail a list of likelihoods for each subject at the observed time points returned by get_likelihood_nofail().
#' @param landmark the prediction landmark.
#' @param horizon the prediction horizon.
#' @param interval the intervals relative to the horizon (e.g. 1/12 if monthly intervals if horizon is in years)
#'
#' @import future.apply
#' @import tidyverse
#'
#' @return A list of length nrow(unique(id)) with likelihood profiles
get_likelihood_profiles_nofail <- function(likelihood_nofail,
                                           landmark,
                                           horizon,
                                           interval) {
  lm <- seq(from = landmark, to = horizon, by = interval)
  obs_l <- future_lapply(likelihood_nofail, function(x) {
    x %>% mutate(likelihoods = exp(cumsum(log(likelihoods))))
  })
  # select the last row and append to build a data.frame of cumulative likelihoods per interval
  Xnew <- future_lapply(obs_l, function(l) {
    temp <- lapply(lm, function(t) {
      temp2 <- l %>% dplyr::filter(time <= t)
      if (nrow(temp2) > 0) {
        temp2 <- temp2 %>%
          dplyr::filter(time == max(time)) %>%
          dplyr::mutate(cutpoint1 = t)
      }
      # else {
      #     temp2 <- data.frame(
      #       id = l$id[1],
      #       time = 0,
      #       likelihoods = 1,
      #       cutpoint1 = t)
      #   }
    })
    do.call(rbind, temp) %>%
      arrange(id) %>%
      dplyr::select(id, cutpoint1, likelihoods)
  })
  Xnew
}


#' Get likelihood profiles for failure
#'
#'Generate a step function of the likelihood at interval times t.
#'
#' @param likelihood_fail a list of likelihoods for each subject at the observed time points
#' @param landmark the prediction landmark time.
#' @param horizon the prediction horizon.
#' @param interval the intervals relative to the horizon (e.g. 1/12 for monthly intervals and horizon is in years)
#'
#' @import future.apply
#' @import tidyverse
#'
#' @return A list of length (horizon - landmark)/interval of lists with nrow(unique(id)) with likelihood profiles
get_likelihood_profiles_fail <- function(likelihood_fail,
                                         landmark,
                                         horizon,
                                         interval) {
  lm <- seq(from = landmark, to = horizon, by = interval)
  # likelihood_fail is a list of lists with length (horizon - landmark + 1)/interval
  # Each of the lists contains a data.frame per patient that needs a cumsum over the likelihood taken.
  list_obs_l <- lapply(likelihood_fail, function(l) {
    future_lapply(seq_along(l), function (x) {
      l[[x]] %>%
        mutate(likelihoods = exp(cumsum(log(likelihoods))))
    })
  })
  # Combine the likelihood patterns for each patient
  obs_l <- lapply(list_obs_l, data.table::rbindlist) %>%
    do.call(rbind, .) %>%
    dplyr::filter(time < pattern) %>%
    split(., by = "id")
  # We now have the patterns per observed time. We need the patterns per possible landmark time (cutpoint 1).
  # dplyr::select the last row and append to build a data.frame of cumulative likelihoods per interval
  Xnew <- future_lapply(obs_l, function(l, ...) {
    temp <- lapply(lm, function(t) {
      temp2 <- l %>%
        dplyr::filter(
          time <= t,
          pattern >= t
        )
      if (nrow(temp2) > 0) {
        temp2 <- temp2 %>%
          dplyr::filter(time == max(time)) %>%
          mutate(cutpoint1 = t)
      }
    })
    do.call(rbind, temp) %>%
      arrange(id) %>%
      dplyr::select(id, cutpoint1, likelihoods, cutpoint2 = pattern)
  })
  Xnew
}

#' Get priors
#'
#' Call flexsurvreg to fit an unconditional survival model and return obtain priors
#'
#' @param data a data.frame with the data for which to predict in the long format (input should be the same as likelihood_profiles).
#' @param time_failure a character with the failure time variable name.
#' @param failure a character with the failure indicator variable name.
#' @param horizon a 1L int with the prediction horizon.
#' @param interval the intervals relative to the horizon (e.g. 1/12 for monthly intervals and horizon is in years)
#' @param dist distribution parameter for the survreg()
#'
#' @import flexsurv
#' @import tidyverse
#' @import survival
#'
#' @return a priori survival and failure probabilities (1-survival) for every interval between min landmark and max horizon
#'
#' @export
get_priors <- function (data, time_failure, failure, horizon, interval, dist = "weibull") {
  surv <- survival::Surv(get(time_failure), get(failure))
  sfit_wb <- flexsurvreg(as.formula(paste(surv,"~ 1")), data = data, dist = dist)
  # get the prior failure probabilities from the weibull model
  pct <- seq(0, 1, by = 0.001)
  pred_wb <- predict(sfit_wb,  type = "quantile", p = pct)$.pred[[1]]
  lm <- (0:(horizon*interval^-1+1))*interval
  prior <- sapply(lm, function (i) {
    pred_wb %>%
      dplyr::filter(.pred <= i) %>%
      summarise(prior_fail = max(.quantile))
  })
  prior <- do.call(rbind, prior)
  prior <- data.frame(
    lm = lm,
    prior_fail = prior,
    prior_survival = 1 - prior,
    row.names = NULL
  )
  prior <- prior %>%
    mutate(prior_survival_interval = lead(prior_survival) / prior_survival,
           prior_fail_interval = 1 - prior_survival_interval) %>%
    dplyr::filter(lm <= horizon) %>%
    dplyr::select(lm, prior_survival_interval, prior_fail_interval)
  prior
}

#' Get predictions for non failure
#'
#' Take a cumsum of the posteriors and return the cumulative incidence of the event between interval and horizon
#'
#' @param likelihood_profiles_fail a list with the likelihood profiles for the failure model for each subject
#'
#' @import tidyverse
#'
#' @return a list of dataframes with the predictions for each subject at every interval up to the horizon
get_predictions_nofail <- function(likelihood_profiles_nofail) {
  negloglik_nofail <- lapply(likelihood_profiles_nofail, setNames, c("id", "cutpoint1", "f_nofail"))
}

#' Get Predictions for Failure
#'
#' Take a cumsum of the posteriors and return the cumulative incidence of the event between interval and horizon
#'Implements equation 7.2 (p155) from Steffen's thesis.
#'
#' @param likelihood_profiles_fail a list with the likelihood profiles for the failure model for each subject
#' @param prior a data frame with priors returned by get_priors()
#'
#' @import tidyverse
#'
#' @return a list of dataframes with the predictions for each subject up to the horizon.
get_predictions_fail <- function (likelihood_profiles_fail, prior) {
  # Join the likelihood and priors data, i.e.
  # f_i(y_i^<=k | F_i = k + 0.5) * P(k <= F_i <= k +1)
  predict_fail <- future_lapply(likelihood_profiles_fail, left_join, prior, by =  c("cutpoint1" = "lm"))
  predict_fail <- future_lapply(predict_fail, function(x) {
    x %>%
      mutate(f_fail_interval = likelihoods * (1-prior_survival_interval)) %>%
      arrange(id, cutpoint1, cutpoint2)})
  future_lapply(predict_fail, function(x) {
    x %>%
      group_by(cutpoint1) %>%
      mutate(sumf = cumsum(f_fail_interval),
             cumfprior = cumprod(prior_survival_interval),
             # Sigma_t=k^h [ f_i(y_i^<=k | F_i = k + 0.5) * P(k <= F_i <= k + 1)] / P(lm <= F_i <= h)
             f_fail = cumsum(f_fail_interval) / (1 - cumprod(prior_survival_interval))) %>%
      ungroup() %>%
      dplyr::select(id, cutpoint1, cutpoint2, f_fail)
  })
}

#' Apply Bayes Rule to get posterior predictions
#'
#' Use the posteriors for failures and non-failure to get predictions with Bayes Rule
#'
#' @param f_fail density for the failures at landmarks upto horizon.
#' @param f_nofail density for the non-failures up to horizon.
#' @param prior the priors returned by get_priors.
#'
#' @import tidyverse
#'
#' @returns: a list with posterior predictions at every landmark
get_posteriors <- function (f_fail, f_nofail, prior) {
  # Combine f_fail and f_nofail
  post <- map2(f_fail, f_nofail, right_join, by = c("id", "cutpoint1"))
  # get P(F>hor|F>=lm) from prior
  post <- future_lapply(post, function(x) {
    x %>% left_join(prior, by = c("cutpoint1" = "lm")) %>%
      group_by(cutpoint1) %>%
      mutate(prior_surv = cumprod(prior_survival_interval),
             prior_fail = 1 - prior_surv,
             post_fail = (f_fail * prior_fail) / (f_fail * prior_fail + f_nofail * prior_surv)) %>%
      ungroup() %>%
      dplyr::select(id, landmark = cutpoint1, horizon=cutpoint2, post_fail)
  })
  post
}

#' Predictions from the multivariate mixed model using a disciminant analysis framework.
#'
#' Calculate the posterior predictions from the likelihood profiles and priors. Implements equation 7.1 (p154) from Steffen's thesis.
#'
#' @param data a data.frame with the data for which to predict in the long format (input should be the same as likelihood_profiles)
#' @param outcomes character vector with the outcomes
#' @param id a character string with the subject id
#' @param fixed_formula_nofail formula for the fixed effects of the nofail model
#' @param random_formula_nofail formula for the random effects of the nofail model,
#' @param random_effects_nofail random effect estimates for the no fail model (see MASS::mvnorm)
#' @param parameters_nofail parameter estimates from the non-failures model (see mmm_model)
#' @param fixed_formula_fail formula for the fixed effects of the fail model,
#' @param random_formula_fail formula for the random effects of the fail model,
#' @param random_effects_fail random effect estimates for the fail model (see MASS::mvnorm),
#' @param parameters_fail parameter estimates from the non-failures model (see mmm_model)
#' @param time a character string with the time variable name in data
#' @param failure a character string with the failure variable name in data
#' @param failure_time a character string with the failure time variable in data
#' @param landmark a numeric value indicating the prediction landmark until which data is available for prediction.
#' @param horizon the prediction horizon (maximum follow-up time)
#' @param interval the intervals relative to the horizon (e.g. 1/12 for monthly intervals and horizon is in years)
#' @param prior the priors returned by get_priors().
#'
#' @import future.apply
#' @import tidyverse
#' @import flexsurv
#'
#' @export
#'
#' @return a data.frame with the posterior predictions for each subject at every landmark up to the horizon

mmm_predictions <- function( data,
                             outcomes,
                             fixed_formula_nofail,
                             random_formula_nofail,
                             random_effects_nofail,
                             parameters_nofail,
                             fixed_formula_fail,
                             random_formula_fail,
                             random_effects_fail,
                             parameters_fail,
                             time,
                             failure,
                             failure_time,
                             prior,
                             id,
                             landmark=0,
                             horizon,
                             interval,
                             ...) {
  n1 <- unique(data[id])
  data <- data %>%
    dplyr::filter(all_of(time) <= landmark)
  n2 <- distinct(data[id])
  if (nrow(n1) != nrow(n2)) {
    dropped <- setdiff(n1, n2) %>% unlist
    message(paste("Dropped", length(dropped), "subjects without data prior to landmark"))
  }
  if (nrow(n2) == 0) {
    stop("There are no observations prior to the landmark. Have you set the argmument `landmark=` correcty")
  }
  likelihood_nofail <- get_likelihood_nofail(data = data,
                                             outcomes = outcomes,
                                             fixed_formula_nofail = fixed_formula_nofail,
                                             random_formula_nofail = random_formula_nofail,
                                             time = time,
                                             parameters_nofail = parameters_nofail,
                                             random_effects_nofail = random_effects_nofail,
                                             id = id,
                                             horizon = horizon)
  likelihood_fail <- get_likelihood_fail(data = data,
                                         outcomes = outcomes,
                                         fixed_formula_fail = fixed_formula_fail,
                                         random_formula_fail = random_formula_fail,
                                         time = time,
                                         failure_time = failure_time,
                                         parameters_fail = parameters_fail,
                                         random_effects_fail = random_effects_fail,
                                         id = id,
                                         landmark = landmark,
                                         interval = interval,
                                         horizon = horizon)
  likelihood_profiles_nofail <- get_likelihood_profiles_nofail(likelihood_nofail = likelihood_nofail,
                                                               landmark = landmark,
                                                               horizon = horizon,
                                                               interval = interval)
  likelihood_profiles_fail <- get_likelihood_profiles_fail(likelihood_fail = likelihood_fail,
                                                           landmark = landmark,
                                                           horizon = horizon,
                                                           interval = interval)
  f_nofail <- get_predictions_nofail(likelihood_profiles_nofail = likelihood_profiles_nofail)
  f_fail <- get_predictions_fail(likelihood_profiles_fail = likelihood_profiles_fail,
                                 prior = prior)
  predictions <- get_posteriors(prior = prior,
                                f_fail = f_fail,
                                f_nofail = f_nofail)
  predictions <- lapply(predictions, nest, data = -id)
  predictions <- lapply(predictions, setNames, c(id, "prediction"))
  predictions <- do.call(rbind, predictions)
  data_trajectory <- data %>%
    group_by(all_of(id)) %>%
    dplyr::filter(all_of(time) <= landmark) %>%
    ungroup() %>%
    dplyr::select(all_of(c(id, time, failure, failure_time, outcomes$outcome))) %>%
    nest(data = -id)
  predictions <- data_trajectory %>% left_join(predictions, by = id)
  if (!is.list(predictions)) {
    predictions = list(predictions)
  }
  predictions
}
