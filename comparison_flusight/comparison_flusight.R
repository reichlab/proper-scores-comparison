# First demo of reading in forecasts and evaluating interval scores

source("comparison_flusight/functions_flusight.R")

# the path where the FluSight forecasts are stored
# adapt to your local file system
path_flusight <- "../../CDClogScore/cdc-flusight-ensemble/"

# example:
# get forecasts from model CUBMA, 2016/2017 , 4wk ahead, US
fc <- get_forecasts(model = "CUBMA",
                    path_flusight = path_flusight,
                    target = "1 wk ahead",
                    location = "US National",
                    season = "2016/2017")
# only really works for week-ahead currently

# evaluate median-and interval score for first week
mis(dens = fc$forecast[1, ], observed = fc$truth[1],
    alpha = 0.2, support = fc$support)