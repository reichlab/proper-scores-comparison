# Generate figure illustrating interval scores for KCDE and SARIMA
# Johannes Bracher, May 2020, johannes.bracher@kit.edu

source("comparison_flusight/functions_flusight.R")


# the path where the FluSight forecasts are stored
# adapt to your local file system
path_flusight <- "../../CDClogScore/cdc-flusight-ensemble/"
list.files(path_flusight)

# example:
# get forecasts from one mode, one location, one target
kcde <- get_forecasts(model = "ReichLab_kcde_backfill_none",
                    path_flusight = path_flusight,
                    target = "1 wk ahead",
                    location = "US National",
                    season = "2016/2017")

sarima <- get_forecasts(model = "ReichLab_sarima_seasonal_difference_TRUE",
                      path_flusight = path_flusight,
                      target = "1 wk ahead",
                      location = "US National",
                      season = "2016/2017")

# perform scoring:
scores_kcde_0.2 <- weighted_interval_score_table(kcde, alpha = c(0.2), weights = 1, detailed = TRUE)
scores_sarima_0.2 <- weighted_interval_score_table(sarima, alpha = c(0.2), weights = 1, detailed = TRUE)

scores_kcde_0.2_1 <- weighted_interval_score_table(kcde, alpha = c(0.2, 1), weights = rep(1, 2), detailed = TRUE)
scores_sarima_0.2_1 <- weighted_interval_score_table(sarima, alpha = c(0.2, 1), weights = rep(1, 2), detailed = TRUE)

alpha_detailed <- c(0.02, 0.05, 1:9/10)

scores_kcde_detailed_u <- weighted_interval_score_table(kcde, alpha = alpha_detailed,
                                                        weights = rep(1, length(alpha_detailed)), detailed = TRUE)
scores_sarima_detailed_u <- weighted_interval_score_table(sarima, alpha = alpha_detailed,
                                                          weights = rep(1, length(alpha_detailed)), detailed = TRUE)

scores_kcde_detailed_w <- weighted_interval_score_table(kcde, alpha = alpha_detailed,
                                                        weights = alpha_detailed, detailed = TRUE)
scores_sarima_detailed_w <- weighted_interval_score_table(sarima, alpha = alpha_detailed,
                                                          weights = alpha_detailed, detailed = TRUE)


# # plotting functions:
# plot_wis_0.2_1 <- function(fc, ylim = c(0, max(fc$wis)), ylab = expression(WIS[list(0, 0.2)]), ...){
#   plot(0.5*(fc$penalty.alpha.1 + fc$penalty.alpha.0.2 + fc$width_pi.alpha.0.2),
#        xlab = "week of season", lwd = 2, type = "h",
#        col = "red", ylim = ylim, ylab = ylab, ...)
#   points(0.5*(fc$width_pi.alpha.0.2 + fc$penalty.alpha.1), lwd = 2, type = "h", col = "lightcoral")
#   points(0.5*(fc$width_pi.alpha.0.2), lwd = 2, type = "h", col = "royalblue1")
#   abline(h = 0)
# }

plot_wis <- function(fc, ylim = c(0, max(fc$wis)), ylab = "WIS", ...){
  plot((fc$weighted_interval_score),
       xlab = "week of season", lwd = 2, type = "h",
       col = "red", ylim = ylim, ylab = ylab, ...)
  points((fc$weighted_width_pi + fc$weighted_penalty_l),
         lwd = 2, type = "h", col = "royalblue1")
  points(fc$weighted_penalty_l,
         lwd = 2, type = "h", col = "orange")
  abline(h = 0)
}

# make plot:
pdf("comparison_flusight/plot_kcde_sarima.pdf", width = 8.5, height = 6)

par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 2), las = 1)
layout(matrix(c(rep(1, 3), rep(2, 3), rep(3, 2),
                rep(1, 3), rep(2, 3), rep(3, 2),
                rep(4, 3), rep(5, 3), rep(3, 2),
                rep(4, 3), rep(5, 3), rep(8, 2),
                rep(6, 3), rep(7, 3), rep(8, 2),
                rep(6, 3), rep(7, 3), rep(8, 2)), nrow = 6, byrow = TRUE))


plot_wis(scores_kcde_0.2, ylim = c(0, 6), ylab = expression(IS[0.2]))
legend("top", legend = c("penalty for exceeding upper limit of 80% PI",
                         "width of 80% PI",
                         "penalty for falling below lower limit of 80% PI"),
       col = c("red", "royalblue1", "orange"), lwd = 2, bty = "n", cex = 0.9)
plot_wis(scores_sarima_0.2, ylim = c(0, 6), ylab = expression(IS[0.2]))

plot(c(mean(scores_kcde_0.2$weighted_interval_score), mean(scores_sarima_0.2$weighted_interval_score)), type = "h",
     ylim = c(0, 1.5), xlim = c(0.5, 2.5), axes = FALSE,
     ylab = expression(average~IS[0.2]),
     col = "red", lwd = 2, xlab = "")
axis(1, at = 1:2, labels = c("KCDE", "SARIMA"), cex.axis = 1)
axis(2); box()
lines(c(mean(scores_kcde_0.2$weighted_width_pi + scores_kcde_0.2$penalty_l.alpha.0.2),
        mean(scores_sarima_0.2$weighted_width_pi + scores_sarima_0.2$penalty_l.alpha.0.2)),
      col = "royalblue1", lwd = 2, type = "h")
lines(c(mean(scores_kcde_0.2$penalty_l.alpha.0.2),
        mean(scores_sarima_0.2$penalty_l.alpha.0.2)),
      col = "orange", lwd = 2, type = "h")
# legend("top", legend = c("penalty for\n non-coverage\n of 80% PI", "average width\n of 80% PI"),
#        col = c("red", "lightcoral", "royalblue1"), lwd = 2, bty = "n", cex = 1)


#--------------------------------------
# time series plot
cex.p <- 0.8

add_pi <- function(l, u, col = "royalblue1"){
  polygon(c(seq_along(l), rev(seq_along(u))),
          c(l, rev(u)), col = col, border = col)
}

plot(kcde$truth, pch = 16, xlab = "week of season", ylab = "%wILI", ylim = c(0, 8),
     main = "", cex = cex.p)
add_pi(scores_kcde_0.2_1$l.alpha.0.2, scores_kcde_0.2_1$u.alpha.0.2)
points(kcde$truth, pch = 16, cex = cex.p)
lines(scores_kcde_0.2_1$u.alpha.1, lty = 3)


legend("topleft", legend = c("median", "80% PI", "observed"),
       lty = c(2, 1, NA), pch = c(NA, NA, 15), bty = "n", cex = 0.9,
       col = c("black", "royalblue1", "black"), lwd = c(1, 4, 1))

plot(sarima$truth, pch = 16, xlab = "week of season", ylab = "%wILI", ylim = c(0, 8),
     main = "", cex = cex.p)
add_pi(scores_sarima_0.2_1$l.alpha.0.2, scores_sarima_0.2_1$u.alpha.0.2)
points(sarima$truth, pch = 16, cex = cex.p)
lines(scores_sarima_0.2_1$u.alpha.1, lty = 3)

#------------------------------------
# Plot for WIS with alpha = c(0.02, 0.05, 0.1, ..., 0.9)




plot_wis(scores_kcde_detailed_w, ylim = c(0, 2), expression(WIS))
legend("top", legend = c("weighted penalty for exceedance of upper limit",
                         "weighted width of PIs",
                         "weighted penalty for falling below lower limit"),
       col = c("red", "royalblue1", "orange"), lwd = 2, bty = "n", cex = 0.9)

plot_wis(scores_sarima_detailed_w, ylim = c(0, 2), expression(WIS))

plot(c(mean(scores_kcde_detailed_w$weighted_interval_score),
       mean(scores_sarima_detailed_w$weighted_interval_score)), type = "h",
     ylim = c(0, 0.5), xlim = c(0.5, 2.5), axes = FALSE,
     ylab = expression(average~WIS),
     col = "red", lwd = 2, xlab = "")
axis(1, at = 1:2, labels = c("KCDE", "SARIMA"), cex.axis = 1)
axis(2); box()
lines(c(mean(scores_kcde_detailed_w$weighted_width_pi + scores_kcde_detailed_w$weighted_penalty_l),
        mean(scores_sarima_detailed_w$weighted_width_pi + scores_sarima_detailed_w$weighted_penalty_l)),
      col = "royalblue1", lwd = 2, type = "h")
lines(c(mean(scores_kcde_detailed_w$weighted_penalty_l),
        mean(scores_sarima_detailed_w$weighted_penalty_l)),
      col = "orange", lwd = 2, type = "h")
# legend("top", legend = c("average penalty for\n non-coverage", "average width\n of PIs"),
#        col = c("red", "royalblue1"), lwd = 2, bty = "n", cex = 1)
dev.off()

# numbers for text:

1 - mean(scores_sarima_0.2$penalty_u.alpha + scores_sarima_0.2$penalty_l.alpha> 0)
1 - mean(scores_kcde_0.2$penalty_u.alpha.0.2 + scores_kcde_0.2$penalty_l.alpha.0.2 > 0)


mean(scores_sarima_detailed_w$weighted_interval_score)
mean(scores_kcde_detailed_w$weighted_interval_score)
