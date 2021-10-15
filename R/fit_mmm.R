# Functions to fit the multivariate mixed models

# support link functions for GLMMadaptive::mixed_model() ----

#' Link function for the combination of a binary and continuous marker
#'
#' @import stats
binary.normal <- function (indicator = stop("'indicator must be specified'")) {
  .indicator <- indicator
  env <- new.env(parent = .GlobalEnv)
  assign(".indicator", indicator, envir = env)
  stats <- make.link("identity")
  log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
    ind <- rep(.indicator, times = length(y)/2)
    eta <- as.matrix(eta)
    out <- eta
    # linear mixed model
    sigma <- exp(phis)
    out[ind == 0, ] <- dnorm(y[ind == 0], eta[ind == 0, ], sigma, log = TRUE)
    # mixed logistic model
    probs <- exp(eta[ind == 1, ]) / (1 + exp(eta[ind == 1, ]))
    out[ind == 1,] <- dbinom(y[ind == 1], 1, probs, log = TRUE)
    out
  }
  structure(list(family = "user binary normal", link = stats$name, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, log_dens = log_dens),
            class = "family")
}
#' Link function for the combination of two continuous markers
#'
#' @import stats
normal <- function () {
  env <- new.env(parent = .GlobalEnv)
  stats <- make.link("identity")
  log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
    eta <- as.matrix(eta)
    out <- eta
    # linear mixed models
    sigma <- exp(phis)
    out <- dnorm(y , eta, sigma, log = TRUE)
    out
  }
  # score_eta_fun <- function (y, mu, phis, eta_zi) {
  #   # the derivative of the log density w.r.t. mu
  #   sigma2 <- exp(phis)^2
  #   y_mu <- y - mu
  #   y_mu / sigma2
  # }
  # score_phis_fun <- function (y, mu, phis, eta_zi) {
  #   sigma <- exp(phis)
  #   y_mu <- y - mu
  #   sigma^{-2} / (1 + y_mu / sigma^2) - 1
  # }
  # environment(log_dens) <- environment(score_eta_fun) <- environment(score_phis_fun) <- env
  structure(list(family = "user defined normal" , link = stats$name, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, log_dens = log_dens),
            class = "family")
}

#' Generate fixed formula for stacked data
#'
#' Takes a list of covariates. Each element of the list are the covariates for a pairwise model
#'
#' @param covars a list of character vectors with the covariates names
#'
#' @export
#'
#' @return a list with the fixed formulas for the pairwise model fitting stage.
#'
#' @examples
#' \dontrun{
#' fixed <- make_fixed_formula(covars = c("time", "sex", "time_failure"))
#' }
make_fixed_formula <- function(covars = NULL) {
  if (!is.null(covars)) {
    covars <- c(sapply(covars, function(x) paste(c(x, "Y1"), collapse = "_")),
                sapply(covars, function(x) paste(c(x, "Y2"), collapse = "_")))
    intercepts <- sapply(c("Y1", "Y2"), function(x) paste(c("intercept",x), collapse = "_"))
    x <- sapply(intercepts, function(x) paste(c(x,"X"), collapse = ":"))
    paste("Y ~ -1 +", paste(c(intercepts, x, covars), collapse = " + "))
  } else {
    intercepts <- sapply(c("Y1", "Y2"), function(x) paste(c("intercept",x), collapse = "_"))
    x <- sapply(intercepts, function(x) paste(c(x,"X"), collapse = ":"))
    paste("Y ~ -1 +", paste(c(intercepts, x), collapse = " + "))
  }
}

#' Generate random effects formula for stacked data
#'
#' Takes a list of covariates for the random slopes.
#' Each element of the list are the random slopes for a pairwise model
#'
#' @param id a character value with the name of the vector of subject identifiers.
#' @param covars a list of character vectors with the covariates names
#'
#' @export
#'
#' @return a list with the fixed formulas for the pairwise model fitting stage.
#'
#' @examples
#' # Random intercept at subject level
#' random <- make_random_formula(id = "id")
#'
#' # Add random slopes
#' \dontrun{
#' random <- make_random_formula(id = "id", covars = "time")
#' }
make_random_formula <- function(id, covars = NULL) {
  if (!is.null(covars)) {
    covars <- c(sapply(covars, function(x) paste(c(x, "Y1"), collapse = "_")),
                sapply(covars, function(x) paste(c(x, "Y2"), collapse = "_")))
    intercepts <- sapply(c("Y1", "Y2"), function(x) paste(c("intercept",x), collapse = "_"))
    paste(paste("~ -1 +", paste(c(intercepts, covars), collapse = " + ")), "|", id)

  } else {
    intercepts <- sapply(c("Y1", "Y2"), function(x) paste(c("intercept",x), collapse = "_"))
    paste(paste("~ -1 +", paste(c(intercepts), collapse = " + ")), "|", id)
  }
}

#' Create random start values
#'
#' Take a model info object created by test_input_datatypes and add random start values
#'
#' @param model_families a list with model inf returned by test_input_datatypes
#' @param n_fixed the number of fixed parameters
#' @param n_random the number of random parameters
#'
#' @import Matrix
#'
#' @return a list with model information including:
#' \itemize{
#'  \item{The model family.}
#'  \item{An indicator for the binary outomce in the binary.normal family.}
#'  \item{Starting values for the model fitting stage.}
#' }
#' @export
get_random_start <- function(model_families, n_fixed, n_random) {
  betas <- rnorm(n_fixed, mean = 0, sd = 0.33)
  D <- as.matrix(
    nearPD(
      matrix(runif(n_random^2), ncol = n_random, nrow = n_random)
    )$mat
  )
  phis <- rnorm(2, mean = 0, sd = 0.33)
  model_families$start_values <- lapply(1:length(model_families$families), function(i) {
    if (model_families$families[[i]] == "binomial") {
      list("betas" = betas,
           "D" = D
      )
    } else if (model_families$families[[i]] == "binary.normal") {
      list("betas" = betas,
           "D" = D,
           "phis" = phis[1]
      )
    } else if (model_families $families[[i]] == "normal") {
      list("betas" = betas,
           "D" = D,
           "phis" = phis
      )
    }
  })
  model_families
}


# Fit a mixed_model for all the frames in stacked data and store the starting values
# optimParallel returns error:
# Error in FGgenerator(par = par, fn = fn, gr = gr, ..., lower = lower,  : is.vector(par) is not TRUE
# Therefore use optim.
# The model fitting within a pair is sequential, since optimParallel throws an error.
# This is a limitation when the number of pairs is less than the number of available workers.
# FIX: add a check for the formula inputs. If a '+' sign is forgotten a opaque error message is returned.
# FIX: add a progress bar during model fitting.

#' Fit pairwise mixed models to get start values for the multivariate mixed model
#'
#' Takes the stacked data and the user defined formulas for the fixed and random effects. Leverages the future.apply package to parallelize along pairs.
#' Using a different number of workers than pairs may result in problems during the optimization process.
#'
#' @param stacked_data a list with the stacked data returned by the stack_data() function.
#' @param fixed character string with the fixed model part. This is passed to as.formula().
#' @param random character string with the random model part. This is passed to as.formula().
#' @param pairs a character matrix with the pairs returned by the make_pairs() function.
#' @param model_families a list with family names and indicators returned by the test_input_datatypes() function.
#' @param ... arguments passed to GLMMadaptive::mixed_model()
#'
#' @import future.apply
#' @import GLMMadaptive
#' @import tidyverse
#'
#' @return a list with starting values to determine the subject level contributions to the derivates.

get_mmm_start_values <- function (stacked_data,
                                  fixed,
                                  random,
                                  pairs,
                                  model_families,
                                  user_initial_values = NULL,
                                  nAGQ = 11,
                                  iter_EM = 30,
                                  iter_qN_outer = 10,
                                  tol1=1e-3,
                                  tol2=1e-4,
                                  tol3=1e-8,
                                  ...) {
  prog <- progressr::progressor(along = nrow(pairs))
  out <- future_lapply(1:nrow(pairs), function(i, ...) {
    prog(sprintf("Fitting pairwise model number %g of %g.", i, nrow(pairs)),
         class = "sticky")
    if (model_families$families[i] == "binary.normal") {
      fit <- mixed_model(fixed = as.formula(fixed),
                         random = as.formula(random),
                         data = stacked_data[[i]],
                         family = binary.normal(indicator = model_families$indicators[[i]]),
                         initial_values = user_initial_values[[i]],
                         nAGQ = nAGQ,
                         n_phis = 1,
                         control = list(optimizer = "optim",
                                        optim_method = "BFGS",
                                        numeric_deriv = "fd",
                                        iter_EM = iter_EM,
                                        iter_qN_outer = iter_qN_outer,
                                        tol1=tol1,
                                        tol2=tol2,
                                        tol3=tol3)
      )
      starting_values <- list(betas = fixef(fit), phis = fit$phis, D = fit$D)
    } else if (model_families$families[i] == "normal") {
      fit <- mixed_model(fixed = as.formula(fixed),
                         random = as.formula(random),
                         data = stacked_data[[i]],
                         family = normal(),
                         initial_values = user_initial_values[[i]],
                         nAGQ = nAGQ,
                         n_phis = 2,
                         control = list(optimizer = "optim",
                                        optim_method = "BFGS",
                                        numeric_deriv = "fd",
                                        iter_EM = iter_EM,
                                        iter_qN_outer = iter_qN_outer,
                                        tol1=tol1,
                                        tol2=tol2,
                                        tol3=tol3)
      )
      starting_values <- list(betas = fixef(fit), phis = fit$phis, D = fit$D)
    } else if (model_families$families[i] == "binomial") {
      fit <- mixed_model(fixed = as.formula(fixed),
                         random = as.formula(random),
                         data = stacked_data[[i]],
                         family = binomial(),
                         initial_values = user_initial_values[[i]],
                         nAGQ = nAGQ,
                         control = list(optimizer = "optim",
                                        optim_method = "BFGS",
                                        numeric_deriv = "fd",
                                        iter_EM = iter_EM,
                                        iter_qN_outer = iter_qN_outer,
                                        tol1=tol1,
                                        tol2=tol2,
                                        tol3=tol3)
      )
      starting_values <- list(betas = fixef(fit), D = fit$D)
    } else {
      stop("Model family not found. Has a list with family for every model been provided?") # Fix: this should be evaluated before fitting is started.
    }
    list("fit" = fit, "starting_values" = starting_values)
  })
  out
}

#' Determine the subject level derivatives.
#'
#' Run mixed models without iterating and return subject level derivatives
#'
#' This function uses the start values determined in by the function get_start_values()
#' to determine the subject level contributions to the derivatives. In turn these are used
#' to average the parameters and standard errors for the full multivariate mixed model.
#'
#' @param id The name of the column with subject ids in stacked_data
#' @param stacked_data A list of data.frames returned by the stack_data() function
#' @param fixed A formula (as character string) for the fixed part of the mixed model. Passed to as.formula().
#' @param random A formula (as character string) for the random part of the mixed model. Passed to as.formula().
#' @param pairs A character matrix with pairs returned by the make_pairs function.
#' @param model_families A list with the model families returned by the test_input_datatypes() function.
#' @param start_values A list (length == length(stacked_data)) with start values returned by the get_start_values() function.
#' @param nAGQ The number of knots to use on the Adapative Gaussian Quadrature (passed to GLMMadaptive::mixed_model()).
#'
#' @import future.apply
#' @import GLMMadaptive
#' @import tidyverse
#'
#' @return A list of length nrow(pairs) each with two elements:
#' \itemize{
#' \item{df_hessians}{a list with data.frames storing the subject level Hessians}
#' \item{df_gradients}{a list with data.frames storing the subject level gradients}
#' }

get_mmm_derivatives <- function (stacked_data, id, fixed, random, pairs, model_families, start_values, nAGQ = 11) {
  stopifnot("id should be a character" = is.character(id))
  prog <- progressr::progressor(along = pairs)
  derivatives <- lapply(1:nrow(pairs), function(i) {
    # message("Getting derivatives for pairwise model: ",i, "out of ", nrow(pairs))
    prog(sprintf("Getting derivatives for pairwise model %g out of %g", i, nrow(pairs)),
         class = "sticky")
    subject_out <- future_lapply(1:nrow(unique(stacked_data[[i]][id])), function(j, ...) {
      subject_data <- stacked_data[[i]][stacked_data[[i]][id] == unique(stacked_data[[i]][id])[j,1],]
      # Data copy to work around optim requiring n > 1 subjects
      subject_data_copy <- subject_data
      subject_data_copy[,id] <- paste0(j,"a")
      subject_data <- rbind(subject_data, subject_data_copy)
      if (model_families$families[i] == "binary.normal") {
        fit <- mixed_model(fixed = as.formula(fixed),
                           random = as.formula(random),
                           data = subject_data,
                           family = binary.normal(indicator = model_families$indicators[[i]]),
                           n_phis = 1,
                           nAGQ = nAGQ,
                           initial_values = start_values[[i]]$starting_values,
                           control = list(optimizer = "optim",
                                          optim_method = "BFGS",
                                          numeric_deriv = "fd",
                                          iter_EM = 0,
                                          iter_qN_outer = 0)
        )
        # Store Hessian (undoing the work around)
        h <- cbind("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
                   fit$Hessian/2)
        g <- c("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
               colSums(fit$score_vect_contributions$score.betas/2),
               colSums(fit$score_vect_contributions$score.D/2),
               colSums(fit$score_vect_contributions$score.phis/2))
      } else if (model_families$families[i] == "normal") {
        fit <- mixed_model(fixed = as.formula(fixed),
                           random = as.formula(random),
                           data = subject_data,
                           family = normal(),
                           n_phis = 2,
                           nAGQ = nAGQ,
                           initial_values = start_values[[i]]$starting_values,
                           control = list(optimizer = "optim",
                                          optim_method = "BFGS",
                                          numeric_deriv = "fd",
                                          iter_EM = 0,
                                          iter_qN_outer = 0)
        )
        h <- cbind("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
                   fit$Hessian/2)
        g <- c("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
               colSums(fit$score_vect_contributions$score.betas/2),
               colSums(fit$score_vect_contributions$score.D/2),
               colSums(fit$score_vect_contributions$score.phis/2))
      } else if (model_families$families[i] == "binomial") {
        fit <- mixed_model(fixed = as.formula(fixed),
                           random = as.formula(random),
                           data = subject_data,
                           family = binomial(),
                           nAGQ = nAGQ,
                           initial_values = start_values[[i]]$starting_values,
                           control = list(optimizer = "optim",
                                          optim_method = "BFGS",
                                          numeric_deriv = "fd",
                                          iter_EM = 0,
                                          iter_qN_outer = 0)
        )
        h <- cbind("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
                   fit$Hessian/2)
        g <- c("id" = as.numeric(unique(stacked_data[[i]][id])[j, 1]),
               colSums(fit$score_vect_contributions$score.betas/2),
               colSums(fit$score_vect_contributions$score.D/2))
      }
      names(g) <- colnames(h)
      list("hessians" = h,
           "gradients" = g)
    })
    list(
      "df_hessians" = do.call(rbind, lapply(subject_out, function(l) l$hessians)),
      "df_gradients" = do.call(rbind, lapply(subject_out, function(l) l$gradients))
    )
  })
  derivatives
}

#' Take the parameter estimates, add labels and stack them into a data.frame
#'
#' Extract the parameter estimates from the start_values and label them back to the names of the original data.
#'
#' @param model_families the model_families object returned by test_input_datatypes()
#' @param pairs A character matrix with pairs returned by the make_pairs function.
#' @param start_values The list of start_values (i.e. parameter estimates) returned by the get_start_values() function.
#'
#' @import tidyverse
#'
#' @return A data.frame with the labelled parameter estimates.

stack_parameters <- function (start_values, pairs, model_families) {
  parameters <- lapply(1:nrow(pairs), function(i) {
    names_fixed <- names(start_values[[i]]$starting_values$betas)
    .parameter_label_y1 <- names_fixed[str_detect(names_fixed, "Y1")]
    .parameter_label_y1 <- str_replace(.parameter_label_y1, "intercept", "")
    if (any(str_detect(.parameter_label_y1, ":X")) == TRUE) {
      .parameter_label_y1[str_detect(.parameter_label_y1, "X")] <-
        paste(rev(unlist(str_split(.parameter_label_y1[str_detect(.parameter_label_y1, "X")], ":"))), collapse = "")
    } else {
      .parameter_label_y1[str_detect(.parameter_label_y1, "X")] <-
        paste(unlist(str_split(.parameter_label_y1[str_detect(.parameter_label_y1, "X")], ":")), collapse = "")
    }
    .parameter_label_y1 <- str_replace(.parameter_label_y1, "X", pairs[i, 2])
    .parameter_label_y1 <- str_replace(.parameter_label_y1, "Y1", pairs[i, 1])
    .parameter_label_y1 <- sapply(1:length(.parameter_label_y1), function(x) paste0("b", x - 1, "_", .parameter_label_y1[x]))
    .parameter_value_y1 <- start_values[[i]]$starting_values$betas[str_detect(names_fixed, "Y1")]
    .parameter_label_y2 <- names_fixed[str_detect(names_fixed, "Y2")]
    .parameter_label_y2 <- str_replace(.parameter_label_y2, "intercept", "")
    if (any(str_detect(.parameter_label_y2, pattern = ":X")) == TRUE) {
      .parameter_label_y2[str_detect(.parameter_label_y2, "X")] <-
        paste(rev(unlist(str_split(.parameter_label_y2[str_detect(.parameter_label_y2, "X")], ":"))), collapse = "")
    } else {
      .parameter_label_y2[str_detect(.parameter_label_y2, "X")] <-
        paste(unlist(str_split(.parameter_label_y2[str_detect(.parameter_label_y2, "X")], ":")), collapse = "")
    }
    .parameter_label_y2 <- str_replace(.parameter_label_y2, "X", pairs[i, 1])
    .parameter_label_y2 <- str_replace(.parameter_label_y2, "Y2", pairs[i, 2])
    .parameter_label_y2 <- sapply(1:length(.parameter_label_y2), function(x) paste0("b", x - 1, "_", .parameter_label_y2[x]))
    .parameter_value_y2 <-  start_values[[i]]$starting_values$betas[str_detect(names_fixed, "Y2")]
    fixed <- data.frame(
      parameter_label = c(.parameter_label_y1, .parameter_label_y2),
      parameter_estimate = c(.parameter_value_y1, .parameter_value_y2),
      row.names = NULL
    )
    # Set names and labels for and get values for random effects
    labels_random <- levels(interaction(
      rownames(start_values[[i]]$starting_values$D),
      colnames(start_values[[i]]$starting_values$D),
      sep = "_x_"))
    labels_random <- t(matrix(labels_random,
                              ncol = ncol(start_values[[i]]$starting_values$D),
                              byrow = FALSE))
    diag(labels_random) <- paste0("var_",diag(labels_random))
    labels_random[lower.tri(labels_random)] <- paste0("cov_", labels_random[lower.tri(labels_random)])
    labels_random <- labels_random[lower.tri(labels_random, diag = TRUE)]
    labels_random <- str_replace_all(labels_random, "Y1", pairs[i,1])
    labels_random <- str_replace_all(labels_random, "Y2", pairs[i,2])
    values_random <- start_values[[i]]$starting_values$D[lower.tri(start_values[[i]]$starting_values$D, diag = TRUE)]
    random <- data.frame(
      parameter_label = labels_random,
      parameter_estimate = values_random,
      row.names = NULL
    )
    # Get residual errors
    resid_labels <- c()
    if (model_families$families[i] == "binary") {
      resid_labels <- NULL
    } else if (model_families$families[i] == "binary.normal") {
      if (all(model_families$indicators[[i]] == c(0, 1 ))) {
        resid_labels <- paste0("resid_", pairs[i, 1])
      } else {
        resid_labels <- paste0("resid_", pairs[i, 2])
      }
    } else if (model_families$families[i] == "normal") {
      resid_labels <- c(paste0("resid_", pairs[i, 1]),
                        paste0("resid_", pairs[i, 2]))
    }
    resid <- data.frame(
      parameter_label = resid_labels,
      parameter_estimate = start_values[[i]]$fit$phis,
      row.names = NULL
    )
    out <- bind_rows(fixed, random, resid)
    rownames(out) <- NULL
    out
  })
  do.call(rbind, parameters)
}

#' Stack gradients into a single data.frame
#'
#' Take the gradients returned by the get_mmm_derivatives function and stack these into a single data.frame for further processing
#'
#' @param derivatives list with derivatives returned by the get_mmm_derivatives() function.
#' @param data The data.frame with original data in long format.
#' @param id A character value with the subject identifier in the orginal data.
#' @param parameter_labels The parameter labels returned by the stack_parameters() function.
#' @param pairs A character matrix with pairs returned by the make_pairs function.
#'
#' @return A data.frame with subject level gradients for every longitudinal outcome

stack_gradients <- function (derivatives, data, id, pairs) {
  gradients <- lapply(1:nrow(pairs), function(i) {
    g <- lapply(1:nrow(unique(data[id])), function(j, ...) {
      data.frame(
        "id" = derivatives[[i]]$df_gradients[j, 1],
        "gradient" = derivatives[[i]]$df_gradients[j, -1],
        row.names = NULL
      )
    })
    do.call(rbind, g)
  })
  gradients <- do.call(rbind, gradients)
  gradients[order(gradients$id),]
}

#' Create the helper 'A'  matrix used to average the parameter estimates.
#'
#' Take parameter labels from MMM object and return a p*q matrix where p is unique(parameters) and q is the total number of parameters (including duplicates).
#'
#' @param parameter_labels A data frame with parameter labels returned by the stack_parameters() function.
#'
#' @return A p*q matrix with ones if a parameter was estimated and zeros if it was not estimated in the pairwise step

create_A_matrix <- function (parameter_labels) {
  A <- matrix(0, nrow = length(unique(parameter_labels)), ncol = length(parameter_labels))
  dimnames(A) <- list(
    unique(parameter_labels),
    parameter_labels
  )
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if (rownames(A)[i] == colnames(A)[j]) {
        A[i,j] <- 1
      }
    }
  }
  A
}

#' Get cross products of the gradients
#'
#' Get the cross product of the subject specific gradient vector, sum and take the average
#'
#' @param data The data.frame with original data in long format.
#' @param id A character value with the subject identifier in the orginal data.
#' @param parameter_labels The parameter labels returned by the stack_parameters() function.
#' @param pairs A character matrix with pairs returned by the make_pairs function.
#' @param gradients A data frame with gradients returned by the stack_gradients() function.
#'
#' @import Matrix
#'
#' @return A matrix with the average gradient for all parameters over all subjects.

get_crossproducts <- function (pairs, parameter_labels,  gradients, data, id) {
  cpgrad <- matrix(0, nrow = length(parameter_labels), ncol = length(parameter_labels),
                   dimnames = list(
                     parameter_labels,
                     parameter_labels
                   )
  )
  # The loop sums all the subject level contributions to the gradient.
  for (i in 1:nrow(unique(data[id]))){
    .id <- unlist(unique(data[paste(id)])[i,])
    gsub <- gradients %>%
      dplyr::filter(id == .id) %>%
      dplyr::select(gradient)
    cpgsub <- tcrossprod(as.matrix(gsub))
    cpgrad <- cpgrad + cpgsub
  }
  # Average the gradients
  cpgrad * (1/nrow(unique(data[paste(id)])))
}

#' Average the Hessians
#'
#' ...
#'
#' @param data The data.frame with original data in long format.
#' @param id A character value with the subject identifier in the orginal data.
#' @param pairs A character matrix with pairs returned by the make_pairs function.
#' @param A the A helper matrix returned by the create_A_matrix function.
#' @param derivatives the derivatives returned by get_mmm_derivatives().
#'
#' @import Matrix
#' @import dplyr
#'
#' @return a Sparse matrix of the MMM Hessian
average_hessians <- function(derivatives, pairs, data, id, A) {
  # Outer loop: create block Sparse matrix of the multivariate mixed model Hessian
  for (j in 1:nrow(pairs)) {
    hess <- as.data.frame(derivatives[[j]]$df_hessians)
    hess_pair <- matrix(0, nrow = ncol(derivatives[[j]]$df_hessians)-1,
                        ncol = ncol(derivatives[[j]]$df_hessians)-1)
    # Inner loop: subject level contributions to the pairwise Hessian
    for (i in 1:nrow(unique(data[id]))) {
      .id <- unlist(unique(data[paste(id)])[i,])
      hsub <- hess %>% dplyr::filter(id == .id)
      hsub <- as.matrix(hsub)[,-1]
      hess_pair <- hess_pair + hsub
    }
    if (j == 1) {
      Hessian <- hess_pair
    } else {
      Hessian <- as.matrix(bdiag(Hessian, hess_pair))
    }
  }
  # Average the Hessian
  Hessian <- Hessian * (1/nrow(unique(data[paste(id)])))
  dimnames(Hessian) <- list(
    colnames(A),
    colnames(A)
  )
  Hessian
}

#' Calculate the covariance matrix for the MMM and return it and the robust standard error
#'
#' Use a sandwich estimator to return the standard errors for the full multivariate mixed model.
#'
#' @param hessian The sparse Hessian matrix for the MMM returned by the average_hessians function().
#' @param cpgradient The cross-product and average of the subject level gradients returned by the get_crossproducts() function.
#' @param parameters The parameter estimates returned by the stack_parameters() function.
#' @param A the A helper matrix returned by the create_A_matrix function().
#' @param data The data.frame with original data in long format.
#' @param id A character value with the subject identifier in the orginal data.
#'
#' @import Matrix
#' @import MASS
#'
#' @return A data.frame with parameter estimates and standard errors for the multivariate mixed model.

# Covariance of sqrt(n) * (theta - theta_0) Fieuws 2006 q(4) p453
get_covariance <- function(hessian, cpgradient, parameters, A, data, id) {
  JKJ <- ginv(hessian) %*% cpgradient %*% ginv(hessian)
  dimnames(JKJ) <- list(
    colnames(A),
    colnames(A)
  )
  # Covariance of (theta - theta_0)
  cov <- JKJ * (1/nrow(unique(data[paste(id)])))
  # Average of the parameters
  A <- A * (1/rowSums(A))
  avg_par <- A %*% parameters$parameter_estimate
  covrobust <- A %*% cov %*% t(A)
  serobust <- sqrt(diag(covrobust))
  out <- data.frame(
    "parameter_name" = rownames(A),
    "parameter_estimate" = avg_par,
    "parameter_std_err" = serobust,
    row.names = NULL
  )
  out
}

#' Create the full mixed model output
#'
#' Take the random effect estimates and return vcov and correlation matrices
#'
#' @param estimates A data.frame with estimates returned by the get_covariance() function.
#'
#' @import tidyverse
#' @import Matrix
#'
#' @return A list of 3 elements including:
#' \describe{
#'   \item{estimates}{A data.frame with the parameter estimates and standard errors}
#'   \item{vcov}{The Variance-Covariance matrix for the random effects}
#'   \item{corr}{A correlation matrix for the random effects}
#'}
get_model_output <- function(estimates) {
  # Sort the fixed effects
  fixef <- estimates %>% dplyr::filter(grepl("^b", parameter_name) == TRUE)
  # Sort the random effects
  varest <- estimates %>% dplyr::filter(grepl("^var", parameter_name) == TRUE)
  covest <- estimates %>% dplyr::filter(grepl("^cov", parameter_name) == TRUE)
  # Sort the residuals
  resid <- estimates %>% dplyr::filter(grepl("resid", parameter_name) == TRUE)
  resid$parameter_estimate <- exp(resid$parameter_estimate)
  resid$parameter_std_err <- exp(resid$parameter_std_err)
  # Covariance matrix
  est_D <- matrix(0, nrow = nrow(varest), ncol = nrow(varest))
  dimnames(est_D) <- list(gsub("^.*_x_","", unique(varest$parameter_name)),
                          gsub("^.*_x_","", unique(varest$parameter_name)))
  diag(est_D) <- varest$parameter_estimate
  est_D[lower.tri(est_D)] <- covest$parameter_estimate
  est_D[upper.tri(est_D)] <- t(est_D)[upper.tri(est_D)]
  # check if est_D is positive definite
  if (all(eigen(est_D)$values > 0) == FALSE) {
    # Compute the nearest positive definite matrix
    est_D <- Matrix::nearPD(est_D)$mat
    est_D <- as.matrix(est_D)
    # Replace the variance and covariance estimates with the nearest PD values
    varest$parameter_estimate <- diag(est_D)
    covest$parameter_estimate <- est_D[lower.tri(est_D)]
    # Send a warning to the user.
    message("The Variance-Covariance matrix was not a positive definite.
            Projecting to the closest positive definite.
            Note that your results may be unexpected.")
  }
  est_Dcorr <- stats::cov2cor(est_D)
  # outcome tables
  out <- list("estimates" = rbind(fixef, varest, covest, resid),
              "vcov" = est_D,
              "corr" = est_Dcorr)
  out
}

#' Fit a multivariate mixed model using the pairwise fitting approach
#'
#' Takes the stacked data and the user defined formulas for the fixed and random effects. Leverages the future.apply package to parallelize along pairs.
#' Using a different number of workers than pairs may result in problems during the optimization process. The output of the pairwise model fitting is combined
#' and any parameters that have been estimated multiple times are averaged. The procedure takes the individual level contributions to the Maximum Likelihood in
#' order to estimate the standard errors for the parameter estimates with a robust sandwich estimator.
#'
#' @param data the original data (in wide, unstacked format).
#' @param stacked_data a list with the stacked data returned by the stack_data() function.
#' @param fixed character string with the fixed model part. This is passed to as.formula().
#' @param random character string with the random model part. This is passed to as.formula().
#' @param pairs a character matrix with the pairs returned by the make_pairs() function.
#' @param model_families a list with family names and indicators returned by the test_input_datatypes() function.
#' @param id The name of the column with subject ids in stacked_data.
#' @param iter_EM (default = 100) EM iterations passed to GLMMadaptive::mixed_model().
#' @param iter_qN_outer (default = 10) outer Quasi Newton iterations passed to GLMMadaptive::mixed_model().
#' @param ncores number of workers for multisession parallel plan.
#' @param ... additional arguments passed to GLMMadaptive::mixed_model().
#'
#' @import tidyverse
#' @import GLMMadaptive
#' @import Matrix
#' @import future.apply
#' @import future
#'
#' @return A list of 3 elements including:
#' \itemize{
#'   \item{estimates}{A data.frame with the parameter estimates and standard errors}
#'   \item{vcov}{The Variance-Covariance matrix for the random effects}
#'   \item{corr}{A correlation matrix for the random effects}
#'}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # enable parallel processing on a Windows machine
#' future::plan(multisession, workers=2)
#'
#' # Fit the model
#' mmm_model(fixed = fixed,
#' random = random,
#' id = "id",
#' data = df,
#' stacked_data = df_stacked,
#' pairs = pairs,
#' model_families = model_info,
#' iter_EM = 100,
#' iter_qN_outer = 30,
#' nAGQ = 7)
#'
#' # Reset the parallel plan back to its default.
#' future::plan("default")
#' }
mmm_model <- function (fixed,
                       random,
                       id,
                       data,
                       stacked_data,
                       pairs,
                       model_families,
                       iter_EM = 100,
                       iter_qN_outer = 10,
                       ncores = NULL,
                       ...) {
  if (!id %in% colnames(data)) {
    stop("Could not find `id` in colnames(`data`). Have you specified the correct subject id?")
  }
  message("Fitting pairwise model")
  start_values <- get_mmm_start_values(stacked_data = stacked_data,
                                       fixed = as.formula(fixed),
                                       random = as.formula(random),
                                       pairs = pairs,
                                       model_families = model_families,
                                       iter_EM = iter_EM,
                                       iter_qN_outer = iter_qN_outer,
                                      # ...
  )
  message("Retrieving derivatives")
  derivatives <- get_mmm_derivatives(stacked_data = stacked_data,
                                     id = id,
                                     fixed = as.formula(fixed),
                                     random = as.formula(random),
                                     pairs = pairs,
                                     model_families = model_families,
                                     start_values = start_values)
  message("Compiling model output.")
  parameters <- stack_parameters(start_values = start_values,
                                 pairs = pairs,
                                 model_families = model_families)
  gradients <- stack_gradients(derivatives = derivatives,
                               data = data,
                               id = id,
                               pairs = pairs)
  A <- create_A_matrix(parameter_labels = parameters$parameter_label)
  cpgradient <- get_crossproducts( pairs = pairs,
                                   parameter_labels = parameters$parameter_label,
                                   gradients = gradients,
                                   data = data,
                                   id = id)
  hessian <- average_hessians(derivatives = derivatives,
                              pairs = pairs,
                              data = data,
                              id = id,
                              A = A)
  mmm_estimates <- get_covariance(hessian = hessian,
                                  cpgradient = cpgradient,
                                  parameters = parameters,
                                  A = A,
                                  data = data,
                                  id = id)
  model <- get_model_output(estimates = mmm_estimates)
  model
}
