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
scores_kcde <- score_forecasts(kcde, miS.alpha = 0.2)
scores_sarima <- score_forecasts(sarima, miS.alpha = 0.2)


# plotting function:
plot_mis <- function(fc, ylim = c(0, max(fc$miS)), ylab = "MIS"){
  plot(fc$miS.abs_error + fc$miS.width_pi + fc$miS.penalty,
       xlab = "week of season", lwd = 2, type = "h",
       col = "red", ylim = ylim, ylab = ylab)
  points(fc$miS.width_pi + fc$miS.abs_error, lwd = 2, type = "h", col = "orange")
  points(fc$miS.width_pi, lwd = 2, type = "h", col = "blue")
  abline(h = 0)
}

# make plot:
pdf("plot_kcde_sarima.pdf", width = 8, height = 4)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), las = 1)
layout(matrix(c(rep(1, 3), rep(2, 3), rep(5, 2),
                rep(3, 3), rep(4, 3), rep(5, 2)), nrow = 2, byrow = TRUE))

plot(kcde$truth, pch = 16, xlab = "week of season", ylab = "%wILI", ylim = c(0, 8),
     main = "KCDE")
lines(scores_kcde$miS.m, lty = 2)
lines(scores_kcde$miS.l, lty = 3)
lines(scores_kcde$miS.u, lty = 3)

legend("topleft", legend = c("median", "80% PI", "observed"),
       lty = c(2, 3, NA), pch = c(NA, NA, 15), bty = "n")


plot(sarima$truth, pch = 16, xlab = "week of season", ylab = "%wILI", ylim = c(0, 8),
     main = "SARIMA")
lines(scores_sarima$miS.m, lty = 2)
lines(scores_sarima$miS.l, lty = 3)
lines(scores_sarima$miS.u, lty = 3)

plot_mis(scores_kcde, ylim = c(0, 16), ylab = expression(MIS[0.2]))
plot_mis(scores_sarima, ylim = c(0, 16), ylab = expression(MIS[0.2]))

plot(c(mean(scores_kcde$miS), mean(scores_sarima$miS)), type = "h",
     ylim = c(0, 3), xlim = c(0.5, 2.5), axes = FALSE,
     ylab = expression(average~MIS[0.2]),
     col = "red", lwd = 2, xlab = "")
axis(1, at = 1:2, labels = c("KCDE", "SARIMA"), cex.axis = 1)
axis(2); box()
lines(c(mean(scores_kcde$miS.abs_error + scores_kcde$miS.width_pi),
        mean(scores_sarima$miS.abs_error + scores_sarima$miS.width_pi)),
      col = "orange", lwd = 2, type = "h")
lines(c(mean(scores_kcde$miS.width_pi),
        mean(scores_sarima$miS.width_pi)),
      col = "blue", lwd = 2, type = "h")
legend("top", legend = c("penalty for\n non-coverage", "absolute error", "width of PI"),
       col = c("red", "orange", "blue"), lwd = 2, bty = "n", cex = 1)
dev.off()
