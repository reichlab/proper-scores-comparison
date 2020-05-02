# Comparison of different scoring rules using FluSight data
# Johannes Bracher, May 2020, johannes.bracher@kit.edu

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

# evaluate median-and interval score for one week
averaged_interval_score(dens = fc$forecast[15, ], observed = fc$truth[15],
                    alpha = c(1), support = fc$support, detailed = TRUE)

log_score(dens = fc$forecast[10, ], observed = fc$truth[10],
          support = fc$support)

crps(dens = fc$forecast[10, ], observed = fc$truth[10], support = fc$support)

# evaluate all scores:
head(log_score_table(fc))
head(averaged_interval_score_table(fc, alpha = c(0.2, 1), detailed = TRUE))

# and run the whole thing for all seasons and models, writing out results:
models <- list.dirs(paste0(path_flusight, "model-forecasts/component-models"),
                    full.names = FALSE, recursive = FALSE)


# get target bounds to know which scores to invlude:
target_bounds <- read.csv(paste0(path_flusight, "writing/comparison/data/all-target-bounds.csv"))

# define alpha levels for detailed mean interval score:
alpha_all <- c(0.02, 0.05, seq(from = 0.1, to = 1, by = 0.1))

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

        # evaluate logS:
        logS_temp <- log_score_table(fc_temp)

        # initialize data.frame to store stuff:
        scores_temp <- logS_temp[, !colnames(logS_temp) %in% c("logS.tolerance", "logS.truncate")]

        # evaluate and add MBlogS:
        MBlogS_temp <- log_score_table(fc_temp, tolerance = 5)
        scores_temp$MBlogS <- MBlogS_temp$logS

        # evaluate and add averaged_interval_score for alpha = 0.2
        miS_temp_0.2 <- averaged_interval_score_table(fc_temp, alpha = 0.2, detailed = FALSE)
        scores_temp[, paste0(c("averaged_penalty", "averaged_width_pi", "averaged_interval_score"), "_0.2")] <-
          miS_temp_0.2[, c("averaged_penalty", "averaged_width_pi", "averaged_interval_score")]

        # evaluate and add averaged_interval_score for alpha = 0.2 and 1, i.e. including 2*AE
        miS_temp_0.2_1 <- averaged_interval_score_table(fc_temp, alpha = c(0.2, 1),
                                                      detailed = TRUE)
        scores_temp[, paste0(c("averaged_penalty", "averaged_width_pi", "averaged_interval_score"), "_0.2_1")] <-
          miS_temp_0.2_1[, c("averaged_penalty", "averaged_width_pi", "averaged_interval_score")]

        miS_temp_all <- averaged_interval_score_table(fc_temp, alpha = alpha_all,
                                                  detailed = FALSE)
        scores_temp[, paste0(c("averaged_penalty", "averaged_width_pi", "averaged_interval_score"), "_all")] <-
          miS_temp_all[, c("averaged_penalty", "averaged_width_pi", "averaged_interval_score")]

        # evaluate and add crps
        crps_temp <- crps_table(fc_temp)
        scores_temp$crps <- crps_temp$crps

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

# models <- models[!models %in% c("FluOutlook_Mech", "FluOutlook_MechAug")]
models <- models[!models %in% c("Delphi_Uniform")]

for(season in c("2016/2017")){
  summary_results[[season]] <- list()
  summary_results[[season]]$logS <- summary_results[[season]]$MBlogS <-
    summary_results[[season]]$crps <- summary_results[[season]]$averaged_interval_score_0.2 <-
    summary_results[[season]]$averaged_interval_score_0.2_1 <-
    summary_results[[season]]$averaged_interval_score_all <-
    matrix(NA, nrow = length(models), ncol = 4,
                                            dimnames = list(models, paste(1:4, "wk ahead")))

  for(model in models){
    all_scores <- read.csv(paste0("comparison_flusight/results/scores_",
                                  gsub(pattern = "/", "-", season), "_", model, ".csv"),
                           stringsAsFactors = FALSE)
    summary_results[[season]]$logS[model, ] <- summarize_scores(all_scores, score = "logS")
    summary_results[[season]]$MBlogS[model, ] <- summarize_scores(all_scores, score = "MBlogS")
    summary_results[[season]]$crps[model, ] <- summarize_scores(all_scores, score = "crps")
    summary_results[[season]]$averaged_interval_score_0.2[model, ] <-
      summarize_scores(all_scores, score = "averaged_interval_score_0.2")
    summary_results[[season]]$averaged_interval_score_0.2_1[model, ] <-
      summarize_scores(all_scores, score = "averaged_interval_score_0.2_1")
    summary_results[[season]]$averaged_interval_score_all[model, ] <-
      summarize_scores(all_scores, score = "averaged_interval_score_all")
  }
}



custom_scatter <- function(score1, score2, ...){
  x <- rowMeans(summary_results$`2016/2017`[[score1]])
  y <- rowMeans(summary_results$`2016/2017`[[score2]])
  plot(x, y,
       pch = 16, cex = 0.8, col = rgb(0.2, 0.2, 0.9, alpha = 0.6), ...)
  legend("bottomleft", paste("corr:", round(cor(x,y), 2)), bty = "n")
}

pdf("comparison_flusight/score_comparison.pdf", width = 7, height = 6)
par(mfcol = c(3, 3), las = 1, mar = c(4, 4.5, 0.5, 0.5))

custom_scatter("logS", "averaged_interval_score_0.2",
     xlab = "logS", ylab = expression(IS[0.2]),
     ylim = yl)

custom_scatter("logS", "averaged_interval_score_all",
               xlab = "logS", ylab = expression(AIS[detailed]),
               ylim = yl)

custom_scatter("logS", "crps",
               xlab = "logS", ylab = "CRPS")


custom_scatter("crps", "averaged_interval_score_0.2",
               xlab = "CRPS", ylab = expression(IS[0.2]),
               ylim = yl)

custom_scatter("crps", "averaged_interval_score_all",
               xlab = "CRPS", ylab = expression(AIS[detailed]),
               ylim = yl)


plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = 0:1, ylim = 0:1)

custom_scatter("averaged_interval_score_all", "averaged_interval_score_0.2",
               xlab = expression(AIS[detailed]), ylab = expression(IS[0.2]),
               ylim = yl)

dev.off()

