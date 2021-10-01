# Model performance -------------------------------------------------------

#' Calculate the Area under the Receiver Operating Characteristics Curve
#'
#' Use the timeROC package to calculate the ROC-AUC based on predictions
#' from the mmm_predictions function.
#'
#' @param predictions predictions data from the mmm_predictions function.
#' @param prediction_landmark a numeric vector with the landmark times to select predictions.
#' @param prediction_horizon the prediction horizon *at* which to calculate the ROC-AUC.
#' @param id a character value for the name of the vector of subject identifiers in the predictions data.
#' @param failure_time a character value for the name of the vector with failure times in the predictions data.
#' @param failure a character value for th name of the vector with failure status.
#'
#' @import timeROC
#' @import tidyverse
#'
#' @details Uses the timeROC package to determine the area under the receiver operating characteristics curve.
#'
#' @export
get_auc <- function (
  predictions,
  prediction_landmark,
  prediction_horizon,
  id,
  failure_time,
  failure
) {
  data <- predictions %>%  # This assumes that predictions are always stored in a list of data tables
    unnest(data) %>%
    group_by(get(id)) %>%
    dplyr::select(all_of(c(failure_time, failure))) %>%
    summarize(across(everything(),first)) %>%
    ungroup()
  predict <- predictions %>%
    unnest(prediction) %>%
    group_by(get(id)) %>%
    dplyr::select(-data) %>%
    dplyr::filter(landmark == prediction_landmark,
                  horizon == prediction_horizon) %>%
    ungroup()
  time <- data %>% dplyr::select(all_of(failure_time)) %>% unlist()
  delta <- data %>% dplyr::select(all_of(failure)) %>% unlist()
  marker <- predict$post_fail
  assertthat::assert_that(
    all(
      length(time)==length(delta),
      length(time)==length(marker)
    ),
    msg = paste("The length of the vectors", failure_time, failure,
                "and predictions do not match.")
  )
  roc <- timeROC(T = time,
                 delta = delta,
                 marker = marker,
                 cause = 1,
                 weighting = "marginal",
                 times = prediction_horizon-0.001, # substract a small amount of time to get result
                 iid = TRUE)
  roc <- cbind(roc$Stats, roc$AUC, roc$inference$vect_sd_1, roc$CumulativeIncidence)
  colnames(roc) <- c("failures", "survivors", "censored", "auc", "auc_se", "cuminc")
  roc <- as.data.frame(roc, row.names = FALSE) %>% mutate(landmark = prediction_landmark)
  roc %>%
    dplyr::filter(!is.na(auc)) %>%
    mutate(horizon = prediction_horizon) %>%
    relocate(landmark, horizon, everything())
}

#' Plot predictions from the multivariate mixed model
#'
#' For a given landmark time, plot the probability of failure conditional on survival
#' and longitudinal biomarker evolutions.
#'
#' @param predictions the predictions generated from mmm_predictions()
#' @param outcomes a character vector with the names of the longitudinal outcomes
#' @param id a character value with the name of the subject identifier
#' @param subject a character value with id for subject of interest
#' @param time a scalar value for the time at which the prediction is plotted
#'
#' @import ggplot2
#' @import gridExtra
#' @import tidyverse
#' @import scales
#' @import tidyr
#'
#' @return plots for survival conditional on the longitudinal biomarker trajectory.
#'
#' @export
plot_predictions <- function (predictions,
                              outcomes,
                              id,
                              subject,
                              time) {
  # FIX: Return marker data to the orginal scale
  # check input
  if((length(subject) == 1) == FALSE) {
    stop("Set `subject` to select a single subject")
  }
  # detect binary outcomes
  is_binary <- sapply(outcomes, function(x, ...) {
    all(
      predictions %>%
        dplyr::select(all_of(id), data) %>%
        tidyr::unnest(data) %>%
        distinct(get(x)) %>%
        unlist() %>%
        sort() == c(0, 1)
    )
  })
  markers <- predictions %>%
    dplyr::select(all_of(id), data) %>%
    dplyr::filter(all_of(id) == subject) %>%
    tidyr::unnest(data)
  # Select prediction data
  pred <- predictions %>%
    dplyr::select(all_of(id), prediction) %>%
    dplyr::filter(all_of(id) == subject) %>%
    tidyr::unnest(prediction)
  # }
  # set-up time axis
  if (nrow(markers) == 0) {
    max_x <- 0
  } else {
    max_x <- markers %>%
      dplyr::select(all_of(time)) %>%
      max()
  }
  # Pivot the data for plotting
  markers_trajectory <- markers %>%
    tidyr::pivot_longer(cols = outcomes, names_to ="outcome")
  markers_trajectory <- markers_trajectory %>%
    dplyr::mutate(is_binary = is_binary[markers_trajectory$outcome])
  # IF binary == TRUE
  binary_plot <- ggplot(
    data = markers_trajectory %>%
      dplyr::filter(
        is_binary == TRUE
        ) %>% as.data.frame(),
    aes(x=time, y=value, col=outcome)) +
    geom_point() +
    geom_line(alpha=0.8) +
    scale_x_continuous(name="",
                       limits=c(0, ceiling(max_x)),
                       breaks=scales::breaks_width(0.5)) +
    scale_y_continuous(
      limits=c(0,1),
      name = "longitudinal outcome") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(vars(outcome))
  # If binary == FALSE
  continuous_plot <- ggplot(
    data=markers_trajectory %>%
      filter(
        is_binary == FALSE
      ) %>% as.data.frame(),
    aes(x=time, y=value, col=outcome)) +
    geom_point() +
    geom_line(alpha = 0.8) +
    scale_x_continuous(name = "",
                       limits = c(0, ceiling(max_x)),
                       breaks = scales::breaks_width(0.5)
    ) +
    scale_y_continuous(name = "longitudinal outcome") +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(vars(outcome),
               scales="free")
  # Plot the predictions
  surv_plot <- pred %>%
    dplyr::filter(landmark==min(landmark)) %>%
    ggplot(aes(x=horizon,y=post_fail)) +
    geom_line(alpha = 0.8)  +
    scale_x_continuous(name = "",
                       breaks = scales::breaks_width(1)) +
    scale_y_continuous(name = "Probability of graft failure",
                       limits = c(0, 1),
                       breaks = scales::breaks_width(0.2),
                       position = "right") +
    theme_bw()
  # Combine plots and return
  gridExtra::grid.arrange(
    grobs=list(binary_plot, continuous_plot, surv_plot),
    widths=c(1,1),
    layout_matrix=rbind(
      c(1,3),
      c(2,3)
    )
  )
}

#' Plot the calibration for MMM
#'
#' Create a grobs with calibration plot for a given landmark and horizon.
#'
#' @param predictions predictions data from the mmm_predictions function.
#' @param prediction_landmark a numeric vector with the landmark times to select predictions
#' @param prediction_horizon the prediction horizon *at* which to calculate the ROC-AUC.
#' @param id a character value for the name of the vector of subject identifiers in the predictions data
#' @param failure_time a character value for the name of the vector with failure times in the predictions data
#' @param failure a character value with the name of the vector of failure indicators.
#' @param predictions the predictions returned by the mmm_predictions function
#'
#' @import tidyverse
#' @import ggplot2
#' @import ggthemes
#' @import gridExtra
#' @import scales
#' @import riskRegression
#'
#' @return a list of grobs
#'
#' @export
plot_calibration <- function(
  predictions,
  prediction_landmark,
  prediction_horizon,
  id,
  failure_time,
  failure
) {
  data <- predictions %>%  # This assumes that predictions are always stored in a list of data tables
    unnest(data) %>%
    group_by(get(id)) %>%
    dplyr::select(all_of(c(failure_time, failure))) %>%
    summarize(across(everything(),first)) %>%
    ungroup()
  predict <- predictions %>%
    unnest(prediction) %>%
    group_by(get(id)) %>%
    dplyr::select(-data) %>%
    dplyr::filter(landmark == prediction_landmark,
                  horizon == prediction_horizon) %>%
    ungroup()
  score <- riskRegression::Score(
    list( "p_fail"=predict$post_fail ),
    formula=Surv(stime, failure) ~ 1,
    data=data,
    conf.inf=FALSE,
    times=prediction_horizon,
    metrics=NULL,
    plots="Calibration")
  plt_calibrate <- riskRegression::plotCalibration(
    score,
    cens.method="local",
    round=TRUE,
    rug=TRUE,
    plot=FALSE)
  plt_calibrate$plotFrames$p_fail %>%
    dplyr::mutate(landmark=prediction_landmark)
}
