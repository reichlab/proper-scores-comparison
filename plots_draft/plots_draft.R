# functions to evaluate logS and MIS for neg bin forecasts:
logS_nb <- function(X, mu, size) dnbinom(X, mu = mu, size = size, log = TRUE)

MIS_nb <- function(X, mu, size, alpha = 0.1){
  l <- qnbinom(alpha/2, mu = mu, size = size)
  u <- qnbinom(1 - alpha/2, mu = mu, size = size)
  m <- qnbinom(0.5, mu = mu, size = size)
  score <- 2*abs(X - m) + (u - l) + 2/alpha*pmax(0, X - u) - 2/alpha*pmin(0, X - l)
  return(score)
}

# Specify values for example:
alpha <- 0.2
supp <- 1:300 # support

# "blue" distribution F:
mu_F <- 60
size_F <- 4
col_F <- "cornflowerblue"
col_G <- rgb(1, 0, 0, 0.5)

# "red" distribution G:
mu_G <- 80
size_G <- 10
col_F2 <- "darkblue"
col_G2 <- "darkred"

y <- 180

d <- dnbinom(supp, mu = mu_F, size = size_F)
pi <- qnbinom(c(alpha/2, 1- alpha/2), mu = mu_F, size = size_F)
m <- qnbinom(0.5, mu = mu_F, size = size_F)

ls <- logS_nb(supp, mu = mu_F, size = size_F)
mis <- MIS_nb(supp, mu = mu_F, size = size_F, alpha = alpha)

d_G <- dnbinom(supp, mu = mu_G, size = size_G)
pi_G <- qnbinom(c(alpha/2, 1 - alpha/2), mu = mu_G, size = size_G)
m_G <- qnbinom(0.5, mu = mu_G, size = size_G)

ls_G <- logS_nb(supp, mu = mu_G, size = size_G)
mis_G <- MIS_nb(supp, mu = mu_G, size = size_G, alpha = alpha)


pdf("plots_draft/fig_logS_miS.pdf", width = 11, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 8), las = 1)

plot(NULL, xlim = c(0, 220), ylim = c(-10, -4), xlab = "y", ylab = "logS")
polygon(400*d - 10, col = col_F, border = col_F)
axis(4, at = c(-10, -8, -6, -4), (c(-10, -8, -6, -4) + 10)/400,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
lines(ls, lwd = 2, col = col_F2)

plot(NULL, xlim = c(0, 220), ylim = c(0, 750), xlab = "y", ylab = expression(MIS[0.2]))
polygon(50000*d, col = col_F, border = col_F)
abline(v = pi, lty = 3, col = col_F2)
abline(v = m, lty = 2, col = col_F2)
lines(mis, lwd = 2, col = col_F2)
text(x = c(pi[1], m, pi[2]), y = 750, labels = c("10%", "50%", "90%"), bty = "o")
axis(4, at = c(-0, 250, 500, 750), c(-0, 250, 500, 750)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
dev.off()


pdf("plots_draft/fig_comparison_tails.pdf", width = 11, height = 4)

par(mfrow = c(1, 2), mar = c(4, 4, 1, 8), las = 1)
plot(NULL, xlim = c(0, 220), ylim = c(-10, -4), xlab = "y", ylab = "logS")
polygon(400*d - 10, col = col_F, border = col_F)
polygon(400*d_G - 10, col = col_G, border = col_G)


axis(4, at = c(-10, -8, -6, -4), (c(-10, -8, -6, -4) + 10)/400,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
lines(ls, lwd = 2, col = col_F2)
lines(ls_G, lwd = 2, col = col_G2)


plot(NULL, xlim = c(0, 220), ylim = c(0, 750), xlab = "y", ylab = expression(MIS[0.2]))
polygon(50000*d, col = col_F, border = col_F)
polygon(50000*d_G, col = col_G, border = col_G)


abline(v = pi, lty = 3, col = col_F2)
abline(v = m, lty = 2, col = col_F2)

abline(v = pi_G, lty = 3, col = col_G2)
abline(v = m_G, lty = 2, col = col_G2)

lines(mis, lwd = 2, col = col_F2)
lines(mis_G, lwd = 2, col = col_G2)

text(x = c(pi[1], m, pi[2]), y = 800, labels = c("5%", "50%", "95%"), bty = "o")
axis(4, at = c(-0, 250, 500, 750), c(-0, 250, 500, 750)/50000,
     col = col_F, col.ticks = col_F)
mtext(side = 4, text =  "predictive probability", las = 0, line = 4.5,
      col = col_F)
dev.off()

