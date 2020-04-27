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
              truth = truths$valid.bin_start_incl))
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

# evaluate log score

#' @param

# evaluate the median-and-interval score

#' @param dens a vector of probabilities
#' @param X the observed value
#' @alpha evaluation is based on (1 - alpha)x100% PI
#' @param support the support underlying the probabilities in dens

mis <- function(dens, observed, alpha, support){
  # get quantiles:
  l <- support[min(which(cumsum(dens) >= alpha/2))]
  m <- support[min(which(cumsum(dens) > 0.5))]
  u <- support[min(which(cumsum(dens) > 1 - alpha/2))]

  abs_error <- abs(observed - m)
  width_pi <- (u - l)
  penalty <- 2/alpha*pmax(0, observed - u) - 2/alpha*pmin(0, observed - l)
  score <-  abs_error + width_pi + penalty

  return(list(l = l,
              m = m,
              u = u,
              observed = observed,
              abs_error = abs_error,
              width_pi = width_pi,
              penalty = penalty,
              score = score))
}