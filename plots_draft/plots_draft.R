# functions to evaluate logS and AIS for neg bin forecasts:
logS_nb <- function(X, mu, size) dnbinom(X, mu = mu, size = size, log = TRUE)

# get functions written originally for FluSight comparison:
source("comparison_flusight/functions_flusight.R")
averaged_interval_score_vect <- function(dens, support, vect_observed, alpha){
  sapply(vect_observed, function(dens, support, obs, alpha){
    averaged_interval_score(dens, support, obs, alpha)$averaged_interval_score
  }, dens = dens, support = support, alpha = alpha
  )
}

# Specify values for example:
alpha.a <- 0.2
alpha.b <- c(0.2, 1)
alpha.c <- c(0.02, 0.05, 1:9/10)

supp <- 1:300 # support

# "blue" distribution F:
mu_F <- 60
size_F <- 4
dens_F <- dnbinom(supp, mu = mu_F, size = size_F)
quantiles_F.a <- qnbinom(c(alpha.a/2, 1 - alpha.a/2), mu = mu_F, size = size_F)
quantiles_F.b <- qnbinom(c(alpha.b/2, 1 - alpha.b/2), mu = mu_F, size = size_F)
quantiles_F.c <- qnbinom(c(alpha.c/2, 1 - alpha.c/2), mu = mu_F, size = size_F)

# compute scores:
logS_F <- logS_nb(supp, mu_F, size_F)
AIS_F.a <- averaged_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.a)
AIS_F.b <- averaged_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.b)
AIS_F.c <- averaged_interval_score_vect(dens = dens_F, support = supp, vect_observed = supp, alpha = alpha.c)

col_F <- "darkolivegreen3"
col_F2 <- "darkolivegreen"



# "red" distribution G:
mu_G <- 80
size_G <- 10
dens_G <- dnbinom(supp, mu = mu_G, size = size_G)
quantiles_G.a <- qnbinom(c(alpha.a/2, 1 - alpha.a/2), mu = mu_G, size = size_G)
quantiles_G.b <- qnbinom(c(alpha.b/2, 1 - alpha.b/2), mu = mu_G, size = size_G)
quantiles_G.c <- qnbinom(c(alpha.c/2, 1 - alpha.c/2), mu = mu_G, size = size_G)

col_G <- rgb(0.8, 0.1, 0.1, 0.5)
col_G2 <- "darkred"

# compute scores:
logS_G <- logS_nb(supp, mu_G, size_G)
AIS_G.a <- averaged_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.a)
AIS_G.b <- averaged_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.b)
AIS_G.c <- averaged_interval_score_vect(dens = dens_G, support = supp, vect_observed = supp, alpha = alpha.c)


pdf("plots_draft/fig_logS_AIS.pdf", width = 9, height = 6)
par(mfrow = c(2, 2), mar = c(4, 4.5, 1, 8), las = 1)

plot(NULL, xlim = c(0, 220), ylim = c(-10, -4), xlab = "y", ylab = "logS")
polygon(400*dens_F - 10, col = col_F, border = col_F)
axis(4, at = c(-10, -8, -6, -4), (c(-10, -8, -6, -4) + 10)/400,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
lines(logS_F, lwd = 2, col = col_F2)

plot(NULL, xlim = c(0, 220), ylim = c(0, 750), xlab = "y", ylab = expression(IS[0.2]))
polygon(50000*dens_F, col = col_F, border = col_F)
abline(v = quantiles_F.a, lty = 3, col = col_F2)
lines(AIS_F.a, lwd = 2, col = col_F2)
text(x = quantiles_F.a, y = 750, labels = c("10%", "90%"), bty = "o")
axis(4, at = c(-0, 250, 500, 750), c(-0, 250, 500, 750)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)


plot(NULL, xlim = c(0, 220), ylim = c(0, 750), xlab = "y", ylab = expression(AIS[list(0,0.2)]))
polygon(50000*dens_F, col = col_F, border = col_F)
abline(v = quantiles_F.b, lty = 3, col = col_F2)
lines(AIS_F.b, lwd = 2, col = col_F2)
text(x = quantiles_F.b, y = 750, labels = c("10%", "50%", "90%", ""), bty = "o")
axis(4, at = c(-0, 250, 500, 750), c(-0, 250, 500, 750)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)


plot(NULL, xlim = c(0, 220), ylim = c(0, 750), xlab = "y", ylab = expression(AIS[detailed]))
polygon(50000*dens_F, col = col_F, border = col_F)
abline(v = quantiles_F.c, lty = 3, col = col_F2)
lines(AIS_F.c, lwd = 2, col = col_F2)
axis(4, at = c(0, 250, 500, 750), c(0, 250, 500, 750)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
dev.off()



pdf("plots_draft/fig_comparison_tails.pdf", width = 11, height = 4.5)

y <- 190

par(mfrow = c(1, 2), mar = c(4, 4, 1, 8), las = 1)
plot(NULL, xlim = c(0, 220), ylim = c(-10, -2), xlab = "y", ylab = "logS")
polygon(400*dens_F - 10, col = col_F, border = col_F)
polygon(400*dens_G - 10, col = col_G, border = col_G)


axis(4, at = c(-10, -8, -6, -4, -2), (c(-10, -8, -6, -4, -2) + 10)/400,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
lines(logS_F, lwd = 2, col = col_F2)
lines(logS_G, lwd = 2, col = col_G2)

abline(v = y, lty = 2)

plot(NULL, xlim = c(0, 220), ylim = c(0, 1000), xlab = "y", ylab = expression(AIS[detailed]))

axis(4, at = c(0, 250, 500, 750, 1000), c(0, 250, 500, 750, 1000)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)

polygon(50000*dens_F, col = col_F, border = col_F)
polygon(50000*dens_G, col = col_G, border = col_G)

lines(AIS_F.c, lwd = 2, col = col_F2)
lines(AIS_G.c, lwd = 2, col = col_G2)

abline(v = y, lty = 2)

dev.off()

# numbers for text:
logS_F[supp == y]
logS_G[supp == y]

AIS_F.c[supp == y]
AIS_G.c[supp == y]

sort(quantiles_F.c, decreasing = TRUE)
sort(quantiles_G.c, decreasing = TRUE)

