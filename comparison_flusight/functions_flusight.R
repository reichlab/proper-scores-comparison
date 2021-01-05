# Functions to handle FluSight forecasts
# Johannes Bracher, johannes.bracher@kit.edu, May 2020

#' custom function to get FluSight forecasts for a given location and target
#'
#' @param model the name of the model (name of the subfolder in which forecasts are stored)
#' @param path the directory in which to look for the subfolder given in model
#' @param location the location for which to read in forecasts
#' @param target the target for which to read in forecasts
#' @param season the season for which to extract forecasts
#'
#' @return a list containing a matrix with forecast densities one row per EW, one column per
#' bin), a vector containing the support underlying the aforementioned matrix,
#' info on season, week, year, target and location, a vector of truth values

get_forecasts <- function(model, path_flusight, target, location, season){

  # path where models are actually stored:
  path_forecast_files <- paste0(path_flusight, "/model-forecasts/component-models")

  # read in all file names:
  fls <- list.files(paste0(path_forecast_files, "/", model)) # get fle names
  fls <- fls[grepl(".csv", fls)]

  # overview of years/weeks/seasons associated with file names:
  overview_fls <- overview_files(fls)
  # restrict to files form respective season:
  relevant_fls <- overview_fls[overview_fls$season == season, ]

  # load the forecasts from the first week to get dimensions:
  dat1 <- read.csv(paste0(path_forecast_files, "/", model, "/", relevant_fls$filename[1]))
  colnames(dat1) <- tolower(colnames(dat1))
  dat1 <- dat1[dat1$target == target &
                 dat1$location == location &
                 dat1$type == "Bin", ]

  # now initialize matrix to store forecasts:
  all_forecasts <- matrix(ncol = nrow(dat1), nrow = nrow(relevant_fls))
  rownames(all_forecasts) <- paste0( relevant_fls$year, "-EW", relevant_fls$week)
  colnames(all_forecasts) <- paste0("bin", as.numeric(as.character(dat1$bin_start_incl)))

  for(i in 1:nrow(relevant_fls)){
    dat_i <- read.csv(paste0(path_forecast_files, "/", model, "/", relevant_fls$filename[i]))
    colnames(dat_i) <- tolower(colnames(dat_i))
    dat_i <- dat_i[dat_i$target == target &
                     dat_i$location == location &
                     dat_i$type == "Bin", ]
    all_forecasts[i, ] <- dat_i$value
  }

  # read in truths:
  truths <- read.csv(paste0(path_flusight, "/scores/target-multivals-20172018.csv"))
  colnames(truths) <- tolower(colnames(truths))
  truths <- truths[truths$target == target &
                       truths$location == location &
                       truths$season == season, ]

  if(nrow(truths) != nrow(all_forecasts)){
    warning("Number of provided forecasts and truths differs.")
  }else{
    if(any(!truths$year == relevant_fls$year) |
       any(!truths$calendar.week == relevant_fls$week) |
       any(!truths$season == relevant_fls$seaons)){
      warning("Truth and forecast times do not match.")
    }
  }

  # get support
  support <- as.numeric(as.character(dat1$bin_start_incl))
  if(any(is.na(support))) support <- as.character(dat1$bin_start_incl)

  return(list(forecast = all_forecasts,
              support = support,
              season = season,
              year = relevant_fls$year,
              week = relevant_fls$week,
              target = target,
              location = location,
              truth = truths$valid.bin_start_incl,
              model = model))
}



# helper function to extract the year, week number and season from a file name
#' @param filename the name of the file
#' @return  a named vector contining the file name, week, year and season
time_from_filename <- function(filename){
  week_and_year <- strsplit(filename, split = "-")[[1]][1:2]

  week <- as.numeric(gsub(pattern = "EW", x = week_and_year[1], replace = ""))
  year <- as.numeric(week_and_year[2])
  season <- ifelse(week < 30,
                   paste0(year - 1, "/", year),
                   paste0(year, "/", year + 1))
  return(c(filename = filename,
           week = week,
           year = year,
           season = season))
}

# helper function to apply time_from_filename to a set of filenames
#' @param filenames a vector of file names
#' @return a data.frame containing the file names, weeks, years and seasons
overview_files <- function(filenames){
  overview <- data.frame(t(sapply(filenames, time_from_filename)))
  overview$year <- as.numeric(as.character(overview$year))
  overview$week <- as.numeric(as.character(overview$week))
  # order by time:
  overview <- overview[order(overview$year, overview$week), ]
  rownames(overview) <- NULL
  overview
}

#' helper function: get quantile from density
#'
#' @param dens the density, i.e. a vector containing probabilities (which sum up to 1)
#' @param support the support of the distribution (i.e. dens[1] = Pr(X = support[1]) etc)
#' @param p probabilities for which to return quantiles
quantile_from_dens <- function(dens, support, p){
  sapply(p, function(p) support[min(which(cumsum(dens) >= p))])
}

#' evaluate the absolute error (of the median)
#'
#' @param dens a vector of probabilities
#' @param the support underlying the probabilities in dens
#' @param observed the observed value
#'
#' @return the absolute error
ae <- function(dens, support, observed){
  m <- quantile_from_dens(dens = dens, support = support, p = 0.5)
  return(abs(observed - m))
}

#' evaluate the CRPS:
#'
#' @param dens a vector of probabilities
#' @param observed the observed value
#' @param support the support underlying the probabilities in dens
#' @param observed the observed value
#'
#' @return the CRPS
crps <- function(dens, support, observed){
  steps0 <- diff(support)
  if(max(steps0) - min(steps0) > 0.0001) stop("crps requires equally spaced support vector.")
  step <- steps0[1]

  cdf <- cumsum(dens)
  cdf_truth <- as.numeric(support >= observed)
  crps <- sum((cdf - cdf_truth)^2)*step # 0.1 is the bin width
  return(crps)
}

# evaluate the mean interval score

#' @param dens a vector of probabilities
#' @param observed the observed value
#' @param alpha a vector of probabilities. Evaluation is based on the central (1 - alpha)x100% PIs. 1 needs to be included to include the absolut error
#' @param weights a vector of weights; if nothing is specified the default of alpha/2 is used for alpha < 1 and 1/2 for the absolute error is used
#' @param support the support underlying the probabilities in dens
#' @param detailed Should a more detailed result be returned?
#'
#' @return a named list containing the score, penalty terms and widths of PI. If detailed == TRUE
#' also split up by alpha value
weighted_interval_score <- function(dens, support, observed, alpha = c(0.1, 0.2, 0.5, 1), weights = NULL, detailed = FALSE){

  if(is.null(weights)){
    weights <- alpha/2
    weights[alpha == 1] <- weights[alpha == 1]/2 # weigh down interval corresponding to absolute error as the IS corresponds to 2 times the AE
  }
  
  if(length(weights) != length(alpha)){
    stop("weights and alpha need to be of the same length.")
  }
  
  alpha_half <- alpha/2
  one_m_alpha_half <- 1 - alpha/2
  
  # compute relevant quantiles:
  m <- quantile_from_dens(dens, support, 0.5)
  l <- quantile_from_dens(dens, support, alpha_half)
  u <- quantile_from_dens(dens, support, one_m_alpha_half)
  
  # compute widths of prediction intervals:
  widths_pi <- u - l
  
  # compute penalties:
  penalties_l <- 2/alpha*pmax(0, l - observed)
  penalties_u <- 2/alpha*pmax(0, observed - u)
  
  # compute interval scores at different levels:
  interval_scores <- widths_pi + penalties_l + penalties_u
  
  # name vectors
  names(l) <- names(u) <- names(widths_pi) <- names(penalties_l) <- names(penalties_u) <-
    names(interval_scores) <- names(alpha) <- names(alpha_half) <-
    names(one_m_alpha_half) <- paste0("alpha.", alpha)
  
  # compute combined score:
  numerator <- length(alpha) - 0.5*as.numeric(1 %in% alpha) # count ae only half if contained
  
  weighted_penalty_l <- sum(weights*penalties_l)/numerator
  weighted_penalty_u <- sum(weights*penalties_u)/numerator
  weighted_width_pi <- sum(weights*widths_pi)/numerator
  weighted_interval_score <- sum(weights*interval_scores)/numerator
  
  if(detailed){
    return(list(
      l = l, u = u, observed = observed,
      width_pi = widths_pi, penalty_l = penalties_l, penalty_u = penalties_u,
      interval_score = interval_scores, weights = weights,
      weighted_penalty_l = weighted_penalty_l,
      weighted_penalty_u = weighted_penalty_u,
      weighted_width_pi = weighted_width_pi,
      weighted_interval_score = weighted_interval_score
    ))
  }else{
    list(weighted_penalty_l = weighted_penalty_l,
         weighted_penalty_u = weighted_penalty_u,
         weighted_width_pi = weighted_width_pi,
         weighted_interval_score = weighted_interval_score)
  }
}

# second implementation to check:
weighted_interval_score2 <- function(dens, support, truth){
  is <- function(l, u, truth, alpha){
    is <- (u - l) + 2/alpha*(truth <= l)*(l - truth) + 2/alpha*(truth >= u)*(truth - u)
    return(is)
  }
  
  vals_alpha <- 1 - c(0, 1:9/10, 0.95, 0.98)
  weights <- vals_alpha/2
  weights[1] <- weights[1]/2 # weight down ae as otherwise counted twice
  
  forecast <- quantile_from_dens(dens, support, c(0.01, 0.025, 1:19/20, 0.975, 0.99))
  
  vals_is <- rep(NA, 12)
  for(i in 1:12){
    l <- forecast[12 - i + 1]
    u <- forecast[12 + i - 1]
    vals_is[i] <- is(l, u, truth, vals_alpha[i])
  } 
  return(sum(weights*vals_is)/11.5)
}

# implementation of linear quantile score to check equivalence:
sum_quantile_score <- function(dens, support, observed){
  quantile_score <- function(value, quantile, observed){
    2*((observed <= value) - quantile)*(value - observed)
  }
  
  quantiles <- c(0.01, 0.025, 1:19/20, 0.975, 0.99)
  forecast <- quantile_from_dens(dens, support, quantiles)

  qs <- numeric(length(quantiles))
  for(i in seq_along(quantiles)){
    qs[i] <- quantile_score(forecast[i], quantiles[i], observed)
  }
  return(sum(qs)/23)
}

# evaluate the logarithmic score

#' @param dens a vector of probabilities
#' @param observed the observed value
#' @param support the support underlying the probabilities in dens
#' @param tolerance tolerance if used in multi-bin version
#' (0 for usual log score, 5 for multi-bin version applied to week-ahead targets)
#' @param truncate value at which log scores are truncated
log_score <- function(dens, observed, support, truncate = -10, tolerance = 0){
  inds_correct <- seq(from = which(support == observed) - tolerance,
                      to = which(support == observed) + tolerance, by = 1)
  inds_correct <- inds_correct[inds_correct > 0]
  logS <- log(sum(dens[inds_correct]))
  if(truncate){
    return(max(truncate, logS))
  }else{
    return(logS)
  }
}

#' Wrapper function to evaluate crps for all rows of a forecast data.frame
#'
#' @param fc the forecast list as returned by get_forecasts
#'
#' @return a data.frame containing various information on the forecasts and the crps

crps_table <- function(fc){
  results <- data.frame(model = fc$model,
                        location = fc$location,
                        season = fc$season,
                        year = fc$year,
                        week = fc$week,
                        target = fc$target,
                        truth = fc$truth,
                        crps = NA)

  for(i in 1:nrow(fc$forecast)){
    results$crps[i] <- crps(dens = fc$forecast[i, ], observed = fc$truth[i],
                            support = fc$support)
  }
  return(results)
}


#' Wrapper function to evaluate absolute error for all rows of a forecast data.frame
#'
#' @param fc the forecast list as returned by get_forecasts
#'
#' @return a data.frame containing various information on the forecasts and the ae

ae_table <- function(fc){
  results <- data.frame(model = fc$model,
                        location = fc$location,
                        season = fc$season,
                        year = fc$year,
                        week = fc$week,
                        target = fc$target,
                        truth = fc$truth,
                        ae = NA)

  for(i in 1:nrow(fc$forecast)){
    results$ae[i] <- ae(dens = fc$forecast[i, ], observed = fc$truth[i],
                            support = fc$support)
  }
  return(results)
}

#' Wrapper function to evaluate logS for all rows of a forecast data.frame
#'
#' @param fc the forecast list as returned by get_forecasts
#' @param truncate Should logS and MBlogS be truncated at -10?
#' @param tolerance the tolerance (in bins) for the MBlogS
#'
#' @return a data.frame containing various information on the forecasts and the logS

log_score_table <- function(fc, truncate = -10, tolerance = 0){
  results <- data.frame(model = fc$model,
                        location = fc$location,
                        season = fc$season,
                        year = fc$year,
                        week = fc$week,
                        target = fc$target,
                        truth = fc$truth,
                        logS.tolerance = tolerance,
                        logS.truncate = truncate,
                        logS = NA)

  for(i in 1:nrow(fc$forecast)){
    results$logS[i] <- log_score(dens = fc$forecast[i, ], observed = fc$truth[i],
                                support = fc$support, tolerance = tolerance,
                                truncate = truncate)
  }
  return(results)
}



#' Wrapper function to evaluate weighted_interval score for all rows of a forecast data.frame
#'
#' @param fc the forecast list as returned by get_forecasts
#' @param alpha a vector pf probabilities. The score the average of the interval
#' score of the (1 - alpha)x100% prediction intervals
#'
#' @return a data.frame containing various information on the forecasts and the weighted_and_interval_score
#'
weighted_interval_score_table <- function(fc, alpha = c(0.1, 0.2, 0.5, 1),
                                          weights = rep(1, length(alpha)), detailed = FALSE){
  results <- data.frame(model = fc$model,
                        location = fc$location,
                        season = fc$season,
                        year = fc$year,
                        week = fc$week,
                        target = fc$target,
                        truth = fc$truth)

  # evaluate for first row and initialize columns:
  temp <- weighted_interval_score(dens = fc$forecast[1, ],  support = fc$support,
                              observed = fc$truth[1], alpha = alpha, weights = weights,
                              detailed = detailed)
  results[, names(unlist(temp))] <- NA
  results[1 , names(unlist(temp))] <- unlist(temp)


  # evaluate remaining rows:
  for(i in 2:nrow(fc$forecast)){
    temp <- weighted_interval_score(dens = fc$forecast[i, ],  support = fc$support,
                                observed = fc$truth[i], alpha = alpha, weights = weights,
                                detailed = detailed)
    results[i, names(unlist(temp))] <- unlist(temp)
  }
  return(results)
}

#' Determine which weeks to keep in the evaluation for a given season, location and target
#'
#' @param target_bounds a table containing the arget bounds, read in from all-target-bounds.csv
#' @param season the season
#' @param location the location
#' @param target the target
#'
#' @return a vector of EWs to be included for evaluation
#'
extract_weeks_evaluated <- function(target_bounds, season, location, target){
  relevant_row <- target_bounds[target_bounds$Season == season &
                                  target_bounds$Location == location &
                                  target_bounds$Target == target, ]
  return(c(relevant_row$start_week:53, # keep 53 in case there is a 53rd week
           1:relevant_row$end_week))
}

#' Compute an averge score from return object of score_forecast
#'
#' @param tab_scores a table containing scores, as returned by score_forecast
#' @param score the name of the score for wich to obtain a summary
#'
#' @return  vector of average scores per target
summarize_scores <- function(tab_scores, score){
  aggr <- aggregate(x = tab_scores[, score], by = list(tab_scores$target), FUN = mean)
  ret <- aggr[, 2]; names(ret) <- aggr[, 1]
  return(ret)
}

#' Compute ransdomized PIT value for one forecast
#'
#' @param dens the discrete density
#' @param support the associated support
#' @param observed the observed value
pit <- function(dens, support, observed){
  sum(dens[support < observed]) + runif(1, 0, 1)*dens[support == observed]
}

#' Compute PIT values for forecast list
#'
#' @param fc the forecast list as returned by get_forecasts
#' @param repetitions How many times should the randomized PIT values be computed?
#'
#' The "splitting" between different bins in case the probability mass assigned to an observed value
#' spans two or more bins is done by repeating the randomized procedure many times rather than
#' computing the values analytically (this was faster to implement)
pit_table <- function(fc, repetitions = 1000){
  pit_vals <- NULL
  for(i in 1:repetitions){
    temp <- numeric(nrow(fc$forecast))
    for(j in 1:nrow(fc$forecast)){
      temp[j] <- pit(dens = fc$forecast[j, ],
                     support = fc$support,
                     observed = fc$truth[j])
    }
    pit_vals <- c(pit_vals, temp)
  }
  return(pit_vals)
}