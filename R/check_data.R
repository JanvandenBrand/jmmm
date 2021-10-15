# Convenience functions to check and wrangle data

#' Make pairs
#'
#' Take a character with the names of the longitudinal outcome variables and make all unique pairs
#'
#' @param outcomes a character vector with the names of longitudinal outcomes
#'
#' @export
#'
#' @returns a matrix of length(outomes)*(length(outcomes)-1)/2 rows and 2 columns with all unique pairs.
#'
#' @examples
#' pairs <- make_pairs(outcomes = c("y1", "y2", "y3", "y4"))
make_pairs <- function(outcomes) {
  N <- length(outcomes)
  Npair <- N * (N-1)/2
  pairs <- matrix(data = 0, nrow = Npair, ncol = 2)
  k <- 1
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      pairs[k,1] <- outcomes[i]
      pairs[k,2] <- outcomes[j]
      k <- k + 1
    }
  }
  pairs
}


#' Test input data types and return model families and indicator for binary outcomes
#'
#' Takes a data.frame with longitudinal data in long format and checks if the longitudinal data are binary or numeric.
#'
#' @param data A data.frame in long format
#' @param pairs A character matrix constructed by make_pairs function.
#'
#' @export
#'
#' @return A list of length nrow(pairs) that includes
#' \describe{
#'   \item{families}{character vector with the custom model family.}
#'   \item{indicator}{a character vector with the indicator for the binary outcome if family = binary.normal}
#' }
#'
#' @examples
#' \dontrun{
#' model_info <- test_input_datatypes(data = df, pairs = pairs)
#' }
test_input_datatypes <- function (data, pairs) {

  # initialize return
  .families <- vector("character", length = nrow(pairs))
  .indicators <- vector("list", length = nrow(pairs))
  # Function body
  for (i in 1:nrow(pairs)) {
    # Assert that y1 is dichotomous
    if (all(sort(unlist(unique(data[paste(pairs[i, 1])]))) == c(0, 1), na.rm = TRUE)) {
      # Assert that y2 is  dichotomous
      if (all(sort(unlist(unique(data[paste(pairs[i, 2])]))) == c(0, 1), na.rm = TRUE)) {
        .families[i] <- "binomial"
        # Assert that y2 is numeric
      } else if (class(unlist(data[paste(pairs[i, 2])])) == "numeric") {
        .families[i] <- "binary.normal"
        .indicators[[i]] <- c(1, 0)
      }
      #  Assert that y1 is numeric
    } else if (class(unlist(data[paste(pairs[i, 1])])) == "numeric") {
      # Assert that y2 is dichotomous
      if (all(sort(unlist(unique(data[paste(pairs[i, 2])]))) == c(0, 1), na.rm = TRUE)) {
        .families[i] <- "binary.normal"
        .indicators[[i]] <- c(0, 1)
        ## Assert that y2 is numeric
      } else if (class(unlist(data[paste(pairs[i, 2])])) == "numeric")  {
        .families[i] <- "normal"
      }
    } else {
      stop("The outcomes variables do not seem to be binary or continuous")
    }
  }
  list("families" = .families, "indicators" = .indicators)
}

#' Return a list with stacked data sets for each of the outcome pairs
#'
#' Takes a data.frame in with longitudinal data to be used in the multivariate mixed model function.
#' Uses the pairwise combination created with the make_pairs() function to unroll the outcome matrix into a vector.
#' The longitudinal predictor and time invariant covariates are processed accordingly.
#'
#' @param data A data.frame with the longitudinal data in long format.
#' @param pairs A character matrix with the pairs, returned from the make_pairs function
#' @param id A character value with the variable name of the subject identifier.
#' @param covars (optional) a character vector with the names of the covariate vectors.
#'
#' @import Matrix
#' @import dplyr
#'
#' @export
#'
#' @return a list of data.frames of length nrow(pairs).
#'
#' @examples
#' \dontrun{
#' df_stacked <- stack_data(
#'   data = df,
#'   id = "id",
#'   pairs = pairs,
#'   covars = c("time", "sex", "time_failure")
#' )
#' }
stack_data <- function (data, id, pairs, covars = NULL) {
  stopifnot("Argument 'id' is missing. Use it to set subject level id." = is.character(id))
  stopifnot("'id' should refer to a character vector" = is.character(unlist(data[id])))
  if (!is.null(covars)) {
    stopifnot("'covars' should be a character vector." = is.character(covars))
  }
  if (!is.null(covars)) {
    # covars should be numeric
    if (all(sapply(data %>% dplyr::select(all_of(covars)), class) == "numeric") == FALSE) {
      stop("stack_data can only handle numeric data: '",
           paste(covars[which(sapply(data %>% dplyr::select(all_of(covars)), class) != "numeric")]),
           "' is not numeric.")
    }
  }
  ids <- rep(unlist(data[id]), each = 2)
  stacked_data <- lapply(1:nrow(pairs), function(i, ...) {
    # i <- 1
    # Unroll Y matrices into vectors
    Y <- data %>%
      dplyr::select(y1 = paste(pairs[i, 1]), y2 = paste(pairs[i, 2]))
    Y <- as.vector(t(as.matrix(Y)))
    Y <- data.frame(id = ids,
                    Y = Y)
    # Create design matrices
    ## The ordering of x1 and x2 is intentional!
    X1 <- data %>%
      dplyr::select(x1 = paste(pairs[i,2]), x2 = paste(pairs[i,1]))
    X1 <- as.vector(t(as.matrix(X1)))
    ## add fixed covariates and time to the design matrices
    X2 <- data %>%
      dplyr::select(all_of(covars))
    X2 <- as.matrix(X2)
    ## add dummies for X
    intercepts <- kronecker(rep(1, times = nrow(data)), diag(1, nrow = 2))
    X2k <- kronecker(X2, diag(1, nrow = 2))
    X2k <- cbind(intercepts, X2k)
    ## join the data
    .stacked_data <- cbind(Y, X1, X2k)
    ## set variable names
    var_names <- c(id, "Y", "X", "intercept_Y1", "intercept_Y2")
    if (!is.null(covars)) {
      .covariate_names <- sapply(covars, function(j) {
        c(paste0(j,"_Y1"), paste0(j,"_Y2"))
      })
      var_names <- c(var_names, paste(.covariate_names))
    }
    names(.stacked_data) <- var_names
    .stacked_data
  })
  stacked_data
}

#' Test to compare stacked outcome data to original outcome data
#'
#' ...
#'
#' @param data data.frame with original data in long format
#' @param stacked_data the list of data.frames returned by stack_data()
#' @param pairs A character matrix with the pairs, returned from the make_pairs function
#'
#' @export
#'
#' @import dplyr
#'
#' @return Prints a list with assertions. All should return TRUE. Otherwise something went wrong.
#'
#' @examples
#' \dontrun{
#' model_info <- test_input_datatypes(data = df, pairs = pairs, stacked_data=df_stacked)
#' }
test_compare_stacked_to_original_data <- function (data, stacked_data, pairs) {
  for (i in 1:nrow(pairs)) {
    current <- unlist(data[paste(pairs[i,1])])
    target <- stacked_data[[i]] %>%
      dplyr::filter(intercept_Y1 == 1) %>%
      dplyr::select(Y)
    check <- cbind(current, target)
    print(paste('Assert: outcome vector matches the original data for', pairs[i, 1], '=',
                all.equal(check[,1], check[,2], na.rm = TRUE, check.attributes = FALSE, use.names = FALSE))
    )
    current <- unlist(data[paste(pairs[i,2])])
    target <- stacked_data[[i]] %>%
      dplyr::filter(intercept_Y1 == 1) %>%
      dplyr::select(X)
    check <- cbind(current, target)
    print(paste('Assert: outcome vector matches the original data for', pairs[i, 2], '=',
                all.equal(check[,1], check[,2], na.rm = TRUE, check.attributes = FALSE, use.names = FALSE))
    )
  }
}

