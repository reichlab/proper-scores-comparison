# First demo of reading in forecasts and evaluating interval scores

source("comparison_flusight/functions_flusight.R")
options(warn = 1)

# the path where the FluSight forecasts are stored
# adapt to your local file system
path_flusight <- "../../CDClogScore/cdc-flusight-ensemble/"

# example:
# get forecasts from one mode, one location, one target
fc <- get_forecasts(model = "CU_RHF_SIRS",
                    path_flusight = path_flusight,
                    target = "1 wk ahead",
                    location = "US National",
                    season = "2016/2017")
# only really works for week-ahead currently

# evaluate median-and interval score for first week
miScore(dens = fc$forecast[1, ], observed = fc$truth[1],
        alpha = 0.2, support = fc$support)

logScore(dens = fc$forecast[1, ], observed = fc$truth[1],
         support = fc$support)

logScore(dens = fc$forecast[1, ], observed = fc$truth[1],
         support = fc$support, tolerance = 1)

# get target bounds to know which scores to invlude:
target_bounds <- read.csv(paste0(path_flusight, "writing/comparison/data/all-target-bounds.csv"))


# evaluate all scores:
head(score_forecasts(fc, miS.alpha = 0.2))
rm(fc)

# and run the whole thing for all seasons and models, writing out results:
models <- list.dirs(paste0(path_flusight, "model-forecasts/component-models"),
                    full.names = FALSE, recursive = FALSE)
all_scores <- NULL

for(model in models){
  for(season in c("2016/2017")){
    cat("Starting", model, "--", season, "\n")

    for(target in paste(1:4, "wk ahead")){
      for(location in c("US National", paste("HHS Region", 1:10))){
        fc_temp <- get_forecasts(model = model,
                                 path_flusight = path_flusight,
                                 target = target,
                                 location = location,
                                 season = season)

        scores_temp <- score_forecasts(fc_temp, miS.alpha = 0.2)

        # restrict to the weeks which were actually evaluated:
        weeks_evaluated <- extract_weeks_evaluated(target_bounds = target_bounds,
                                                   season = season,
                                                   location = location,
                                                   target = target)
        scores_temp <- subset(scores_temp, week %in% weeks_evaluated)

        if(is.null(all_scores)){
          all_scores <- scores_temp
        }else{
          all_scores <- rbind(all_scores, scores_temp)
        }
      }
    }

    write.csv(all_scores, file = paste0("comparison_flusight/results/scores_",
                                        gsub(pattern = "/", "-", season), "_", model, ".csv"),
              row.names = FALSE)
    all_scores <- NULL

  }
}


# analyse the scores:
summary_results <- list()

models <- models[!models %in% c("FluOutlook_Mech", "FluOutlook_MechAug")]

for(season in c("2016/2017")){
  summary_results[[season]] <- list()
  summary_results[[season]]$logS <- summary_results[[season]]$MBlogS <-
    summary_results[[season]]$miS <- matrix(NA, nrow = length(models), ncol = 4,
                                            dimnames = list(models, paste(1:4, "wk ahead")))

  for(model in models){
    all_scores <- read.csv(paste0("comparison_flusight/results/scores_",
                                  gsub(pattern = "/", "-", season), "_", model, ".csv"),
                           stringsAsFactors = FALSE)
    summary_results[[season]]$logS[model, ] <- summarize_scores(all_scores, score = "logS")
    summary_results[[season]]$MBlogS[model, ] <- summarize_scores(all_scores, score = "MBlogS")
    summary_results[[season]]$miS[model, ] <- summarize_scores(all_scores, score = "miS")
  }
}

par(mfrow = 1:2)
plot(summary_results$`2016/2017`$logS[, 1], summary_results$`2016/2017`$miS[, 1],
     main = "1 week ahead", xlab = "logS", ylab = "miS")

plot(summary_results$`2016/2017`$MBlogS[, 1], summary_results$`2016/2017`$miS[, 1],
     main = "1 week ahead", xlab = "MBlogS", ylab = "miS")


