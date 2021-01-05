# Comparison of different scoring rules using FluSight data
# Johannes Bracher, May 2020, johannes.bracher@kit.edu

source("comparison_flusight/functions_flusight.R")
options(warn = 1)

# the path where the FluSight forecasts are stored
# adapt to your local file system
path_flusight <- "../cdc-flusight-ensemble/"

# example:
# get forecasts from one mode, one location, one target
fc <- get_forecasts(model = "CU_RHF_SIRS",
                    path_flusight = path_flusight,
                    target = "1 wk ahead",
                    location = "US National",
                    season = "2016/2017")
# only really works for week-ahead currently

# evaluate weighted interval score for one week
ind <- 16

weighted_interval_score(dens = fc$forecast[ind, ], observed = fc$truth[ind],
                        alpha = c(0.02, 0.05, 1:10/10), support = fc$support, detailed = TRUE)

# check against two simpler/less detailed implementations:
weighted_interval_score2(dens = fc$forecast[ind, ], support = fc$support, truth = fc$truth[ind])
sum_quantile_score(dens = fc$forecast[ind, ], support = fc$support, observed = fc$truth[ind])

log_score(dens = fc$forecast[ind, ], observed = fc$truth[ind],
          support = fc$support)

crps(dens = fc$forecast[ind, ], observed = fc$truth[ind], support = fc$support)
# CRPS is approximated to WIS with a low of intervals:
weighted_interval_score(dens = fc$forecast[ind, ], observed = fc$truth[ind],
                        alpha = 1:1000/1000, support = fc$support, detailed = FALSE)


ae(dens = fc$forecast[ind, ], observed = fc$truth[ind], support = fc$support)


# evaluate all scores:
head(log_score_table(fc))
head(weighted_interval_score_table(fc, alpha = c(0.2, 1), detailed = TRUE))

# and run the whole thing for all seasons and models, writing out results:
models <- list.dirs(paste0(path_flusight, "model-forecasts/component-models"),
                    full.names = FALSE, recursive = FALSE)


# get target bounds to know which scores to invlude:
target_bounds <- read.csv(paste0(path_flusight, "writing/comparison/data/all-target-bounds.csv"))

# define alpha levels for detailed mean interval score:
alpha_detailed <- c(0.02, 0.05, seq(from = 0.1, to = 1, by = 0.1))

# create directory to store results:
dir.create("comparison_flusight/results")

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

        # evaluate and add weighted_interval_score for alpha = 0.2
        miS_temp_0.2 <- weighted_interval_score_table(fc_temp, alpha = 0.2, detailed = FALSE, weights = 1)
        scores_temp[, paste0(c("weighted_penalty_l", "weighted_penalty_u",
                               "weighted_width_pi", "weighted_interval_score"), "_0.2")] <-
          miS_temp_0.2[, c("weighted_penalty_l", "weighted_penalty_u",
                           "weighted_width_pi", "weighted_interval_score")]

        # # evaluate and add weighted_interval_score for alpha = 0.2 and 1, i.e. including 2*AE
        # miS_temp_0.2_1 <- weighted_interval_score_table(fc_temp, alpha = c(0.2, 1),
        #                                               detailed = TRUE)
        # scores_temp[, paste0(c("weighted_penalty", "weighted_width_pi", "weighted_interval_score"), "_0.2_1")] <-
        #   miS_temp_0.2_1[, c("weighted_penalty", "weighted_width_pi", "weighted_interval_score")]

        # evaluate and add unweighted detailed score:
        miS_temp_detailed_u <- weighted_interval_score_table(fc_temp, alpha = alpha_detailed,
                                                             weights = rep(1, length(alpha_detailed)), detailed = FALSE)
        scores_temp[, paste0(c("weighted_penalty_l", "weighted_penalty_u",
                               "weighted_width_pi", "weighted_interval_score"), "_detailed_u")] <-
          miS_temp_detailed_u[, c("weighted_penalty_l", "weighted_penalty_u",
                                  "weighted_width_pi", "weighted_interval_score")]

        # evaluate detailed score with alpha as weights:
        miS_temp_detailed_w <- weighted_interval_score_table(fc_temp, alpha = alpha_detailed,
                                                             weights = alpha_detailed/2, detailed = FALSE)
        scores_temp[, paste0(c("weighted_penalty_l", "weighted_penalty_u",
                               "weighted_width_pi", "weighted_interval_score"), "_detailed_w")] <-
          miS_temp_detailed_w[, c("weighted_penalty_l", "weighted_penalty_u",
                                  "weighted_width_pi", "weighted_interval_score")]

        # evaluate and add crps
        crps_temp <- crps_table(fc_temp)
        scores_temp$crps <- crps_temp$crps

        # evaluate and add absolute error:
        ae_temp <- ae_table(fc_temp)
        scores_temp$ae <- ae_temp$ae

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
    summary_results[[season]]$neg_logS <- summary_results[[season]]$neg_MBlogS <-
    summary_results[[season]]$crps <- summary_results[[season]]$weighted_interval_score_0.2 <-
    summary_results[[season]]$weighted_interval_score_0.2_1 <-
    summary_results[[season]]$weighted_interval_score_detailed_u <-
    summary_results[[season]]$weighted_interval_score_detailed_w <-
    summary_results[[season]]$ae <-
    matrix(NA, nrow = length(models), ncol = 4,
           dimnames = list(models, paste(1:4, "wk ahead")))

  for(model in models){
    all_scores <- read.csv(paste0("comparison_flusight/results/scores_",
                                  gsub(pattern = "/", "-", season), "_", model, ".csv"),
                           stringsAsFactors = FALSE)
    summary_results[[season]]$logS[model, ] <- summarize_scores(all_scores, score = "logS")
    summary_results[[season]]$neg_logS[model, ] <- -1*summarize_scores(all_scores, score = "logS")

    summary_results[[season]]$MBlogS[model, ] <- summarize_scores(all_scores, score = "MBlogS")
    summary_results[[season]]$neg_MBlogS[model, ] <- -1*summarize_scores(all_scores, score = "MBlogS")

    summary_results[[season]]$crps[model, ] <- summarize_scores(all_scores, score = "crps")
    summary_results[[season]]$weighted_interval_score_0.2[model, ] <-
      summarize_scores(all_scores, score = "weighted_interval_score_0.2")
    # summary_results[[season]]$weighted_interval_score_0.2_1[model, ] <-
    #   summarize_scores(all_scores, score = "weighted_interval_score_0.2_1")
    summary_results[[season]]$weighted_interval_score_detailed_u[model, ] <-
      summarize_scores(all_scores, score = "weighted_interval_score_detailed_u")
    summary_results[[season]]$weighted_interval_score_detailed_w[model, ] <-
      summarize_scores(all_scores, score = "weighted_interval_score_detailed_w")
    summary_results[[season]]$ae[model, ] <-
      summarize_scores(all_scores, score = "ae")
  }
}


# custom function for scatter plot:
custom_scatter <- function(score1, score2, show_diag = FALSE,...){
  x <- rowMeans(summary_results$`2016/2017`[[score1]])
  y <- rowMeans(summary_results$`2016/2017`[[score2]])

  cols <- ifelse(grepl("FluOutlook_Mech", names(x)),
                 rgb(0.9, 0.6, 0, alpha = 0.85),
                 rgb(0.2, 0.2, 0.9, alpha = 0.85))

  plot(x, y, col = NA, ...)
  if(show_diag) abline(0:1, col = "lightgrey")
  points(x, y, col = cols, pch = 16, cex = 0.8)
  legend("bottomright", paste("corr:", round(cor(x,y), 2)), bty = "n")
}

yl_is <- c(0, 8)
yl_wis <- c(0, 1)

pdf("comparison_flusight/score_comparison.pdf", width = 7, height = 10)
par(mfrow = c(5, 2), las = 1, mar = c(4, 4.5, 1.5, 0.5))

# row for log score:
custom_scatter("neg_logS", "weighted_interval_score_0.2",
               xlab = "-logS", ylab = expression(IS[0.2]),
               ylim = yl_is)

custom_scatter("neg_logS", "weighted_interval_score_detailed_w",
               xlab = "-logS", ylab = expression(WIS),
               ylim = yl_wis)


# row for MBlogS:
custom_scatter("neg_MBlogS", "weighted_interval_score_0.2",
               xlab = "-MBlogS", ylab = expression(IS[0.2]),
               ylim = yl_is)

custom_scatter("neg_MBlogS", "weighted_interval_score_detailed_w",
               xlab = "-MBlogS", ylab = expression(WIS),
               ylim = yl_wis)


# row for CRPS:
custom_scatter("crps", "weighted_interval_score_0.2",
               xlab = "CRPS", ylab = expression(IS[0.2]),
               ylim = yl_is, xlim = yl_wis)

custom_scatter("crps", "weighted_interval_score_detailed_w",
               xlab = "CRPS", ylab = expression(WIS),
               ylim = yl_wis, xlim = yl_wis, show_diag = TRUE)

# row for AE:
custom_scatter("ae", "weighted_interval_score_0.2",
               xlab = "AE", ylab = expression(IS[0.2]),
               ylim = yl_is, xlim = yl_wis)

custom_scatter("ae", "weighted_interval_score_detailed_w",
               xlab = "AE", ylab = expression(WIS),
               ylim = yl_wis, xlim = yl_wis, show_diag = TRUE)

# row for IS:
plot(NULL, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")

custom_scatter("weighted_interval_score_0.2", "weighted_interval_score_detailed_w",
               xlab = expression(IS[0.2]), ylab = expression(WIS),
               ylim = yl_wis, xlim = yl_is)


dev.off()

