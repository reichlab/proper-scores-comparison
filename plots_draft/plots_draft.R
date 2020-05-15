# get functions written originally for FluSight comparison:
source("comparison_flusight/functions_flusight.R")

# functions to evaluate logS, CRPS and WIS for neg bin forecasts:
logS_nb <- function(X, mu, size) dnbinom(X, mu = mu, size = size, log = TRUE)

weighted_interval_score_vect <- function(dens, support, vect_observed, alpha, weights = alpha/2){
  sapply(vect_observed, function(dens, support, obs, alpha, weights){
    weighted_interval_score(dens = dens, support = support,
                            obs = obs, alpha = alpha,
                            weights = weights)$weighted_interval_score
  }, dens = dens, support = support, alpha = alpha, weights = weights
  )
}

crps_vect <- function(dens, support, vect_observed){
  sapply(vect_observed, function(dens, support, obs, alpha){
    crps(dens = dens, support = support, observed = obs)
  }, dens = dens, support = support
  )
}

log_score_vect <- function(dens, support, vect_observed, tolerance = 0, truncate = -10){
  sapply(vect_observed, function(dens, support, obs, alpha, tolerance, truncate){
    log_score(dens = dens, support = support, observed = obs,
              tolerance = tolerance, truncate = truncate)
  }, dens = dens, support = support, tolerance = tolerance, truncate = truncate
  )
}

# Specify values for example:
alpha.a <- 0.2
alpha.b <- c(0.2, 1)
alpha.c <- c(0.1, 0.4, 0.7, 1)
alpha.d <- c(0.02, 0.05, 1:9/10)

supp <- 1:300 # support

# "blue" distribution F:
mu_F <- 60
size_F <- 4
dens_F <- dnbinom(supp, mu = mu_F, size = size_F)
m_F <- quantile_from_dens(dens_F, supp, 0.5)
quantiles_F.a <- qnbinom(c(alpha.a/2, 1 - alpha.a/2), mu = mu_F, size = size_F)
quantiles_F.b <- qnbinom(c(alpha.b/2, 1 - alpha.b/2), mu = mu_F, size = size_F)
quantiles_F.c <- qnbinom(c(alpha.c/2, 1 - alpha.c/2), mu = mu_F, size = size_F)
quantiles_F.d <- qnbinom(c(alpha.c/2, 1 - alpha.d/2), mu = mu_F, size = size_F)


# compute scores:
logS_F <- log_score_vect(dens = dens_F, support = supp, vect_observed = supp, truncate = -100)
MBlogS_F <- log_score_vect(dens = dens_F, support = supp, vect_observed = supp, tolerance = 5)
crps_F <- crps_vect(dens = dens_F, support = supp, vect_observed = supp)
ae_F <- ae(dens = dens_F, support = supp, observed = supp)

WIS_F.a <- weighted_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.a,
                                        weights = 1)
WIS_F.b <- weighted_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.b)
WIS_F.c <- weighted_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.c)
WIS_F.d <- weighted_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.d)


col_F <- "darkolivegreen3"
col_F2 <- "darkolivegreen"



# "red" distribution G:
mu_G <- 80
size_G <- 10
dens_G <- dnbinom(supp, mu = mu_G, size = size_G)
m_G <- quantile_from_dens(dens_G, supp, 0.5)
quantiles_G.a <- qnbinom(c(alpha.a/2, 1 - alpha.a/2), mu = mu_G, size = size_G)
quantiles_G.b <- qnbinom(c(alpha.b/2, 1 - alpha.b/2), mu = mu_G, size = size_G)
quantiles_G.c <- qnbinom(c(alpha.c/2, 1 - alpha.c/2), mu = mu_G, size = size_G)
quantiles_G.d <- qnbinom(c(alpha.d/2, 1 - alpha.c/2), mu = mu_G, size = size_G)


col_G <- rgb(0.8, 0.1, 0.1, 0.5)
col_G2 <- "darkred"

# compute scores:
logS_G <- log_score_vect(dens = dens_G, support = supp, vect_observed = supp, truncate = -100)
MBlogS_G <- log_score_vect(dens = dens_G, support = supp, vect_observed = supp,
                           tolerance = 5, truncate = -100)
crps_G <- crps_vect(dens = dens_G, support = supp, vect_observed = supp)
ae_G <- ae(dens = dens_G, support = supp, observed = supp)


WIS_G.a <- weighted_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.a,
                                        weights = 1)
WIS_G.b <- weighted_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.b)
WIS_G.c <- weighted_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.c)
WIS_G.d <- weighted_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.d)

green_axis <- function(at, labels, cex = 0.7){
  axis(4, at = at, labels = rep("", 4), col = col_F, col.ticks = col_F)
  mtext(labels, side = 4, at = at, line = 1, cex = cex, col = col_F)
  mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
        col = col_F, cex = cex)
}


pdf("plots_draft/fig_logS_WIS.pdf", width = 7.5, height = 6)
par(mfrow = c(3, 2), mar = c(4, 4.5, 1, 8), las = 1)
yl <- c(0, 750)

plot(NULL, xlim = c(0, 220), ylim = c(-10, -4), xlab = "y", ylab = "logS")
polygon(400*dens_F - 10, col = col_F, border = col_F)
green_axis(at = c(-10, -8, -6, -4), labels = (c(-10, -8, -6, -4) + 10)/400)
lines(logS_F, lwd = 2, col = col_F2)


plot(NULL, xlim = c(0, 220), ylim = c(0, 750/4), xlab = "y", ylab = "CRPS")
polygon(50000*dens_F/4, col = col_F, border = col_F)
lines(crps_F, lwd = 2, col = col_F2)
green_axis(at = c(-0, 250, 500, 750)/4, labels = c(-0, 250, 500, 750)/50000)


plot(NULL, xlim = c(0, 220), ylim = c(0, 750/4), xlab = "y", ylab = "absolute error")
polygon(50000*dens_F/4, col = col_F, border = col_F)
lines(ae_F, lwd = 2, col = col_F2)
abline(v = m_F, lty = 3, col = col_F2)
text(x = m_F, y = 750/4 - 10, labels = c("50%"), bty = "o")
green_axis(at = c(-0, 250, 500, 750)/4, labels = c(-0, 250, 500, 750)/50000)



plot(NULL, xlim = c(0, 220), ylim = yl/4, xlab = "y", ylab = expression(WIS^"*"))
polygon(50000*dens_F/4, col = col_F, border = col_F)
abline(v = quantiles_F.c, lty = 3, col = col_F2)
lines(WIS_F.c, lwd = 2, col = col_F2)
text(x = quantiles_F.c[c(1, 3, 7, 5)], y = 170, labels = c("5%", "35%", "65%", "95%"),
     bty = "o", cex = 0.8)
text(x = quantiles_F.c[c(2, 4, 6)], y = 150, labels = c("20%", "50%", "80%"),
     bty = "o", cex = 0.8)

green_axis(at = c(-0, 250, 500, 750)/4, labels = c(-0, 250, 500, 750)/50000)


plot(NULL, xlim = c(0, 220), ylim = yl, xlab = "y", ylab = expression(IS[0.2]))
polygon(50000*dens_F, col = col_F, border = col_F)
abline(v = quantiles_F.a, lty = 3, col = col_F2)
lines(WIS_F.a, lwd = 2, col = col_F2)
text(x = quantiles_F.a, y = 750 - 40, labels = c("10%", "90%"), bty = "o")
green_axis(at = c(-0, 250, 500, 750), labels = c(-0, 250, 500, 750)/50000)


plot(NULL, xlim = c(0, 220), ylim = yl/4, xlab = "y", ylab = expression(WIS))
polygon(50000*dens_F/4, col = col_F, border = col_F)
abline(v = quantiles_F.d, lty = 3, col = col_F2)
lines(WIS_F.d, lwd = 2, col = col_F2)
green_axis(at = c(-0, 250, 500, 750)/4, labels = c(-0, 250, 500, 750)/50000)

dev.off()



pdf("plots_draft/fig_comparison_tails.pdf", width = 10, height = 3.2)

y <- 190

par(mfrow = c(1, 2), mar = c(4, 4, 1, 8), las = 1)
plot(NULL, xlim = c(0, 220), ylim = c(-10, -4), xlab = "y", ylab = "logS")
polygon(400*dens_F - 10, col = col_F, border = col_F)
polygon(400*dens_G - 10, col = col_G, border = col_G)
green_axis(at = c(-10, -8, -6, -4), labels = (c(-10, -8, -6, -4) + 10)/400, cex = 1)

lines(logS_F, lwd = 2, col = col_F2)
lines(logS_G, lwd = 2, col = col_G2)

abline(v = y, lty = 2)

plot(NULL, xlim = c(0, 220), ylim = yl, xlab = "y", ylab = expression(WIS))
green_axis(at = c(-0, 250, 500, 750), labels = c(-0, 250, 500, 750)/50000, cex = 1)


polygon(50000*dens_F, col = col_F, border = col_F)
polygon(50000*dens_G, col = col_G, border = col_G)

lines(WIS_F.d, lwd = 2, col = col_F2)
lines(WIS_G.d, lwd = 2, col = col_G2)

abline(v = y, lty = 2)

dev.off()

# numbers for text:
logS_F[supp == y]
logS_G[supp == y]

WIS_F.c[supp == y]
WIS_G.c[supp == y]

sort(quantiles_F.c, decreasing = TRUE)
sort(quantiles_G.c, decreasing = TRUE)

