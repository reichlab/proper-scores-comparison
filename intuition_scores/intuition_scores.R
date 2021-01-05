# setwd("/home/johannes/Documents/COVID/proper-scores-comparison/intuition_scores")

library(pBrackets)

col_F <- "darkolivegreen3"
col_F2 <- "darkolivegreen"

# define distributions:
mu <- 60
size <- 4
support <- 0:200
y <- 55

dens_F <- dnbinom(support, mu = mu, size = size)
cdf_F <- pnbinom(support, mu = mu, size = size)
q0.1 <- qnbinom(p = 0.1, mu = mu, size = size)
q0.9 <- qnbinom(p = 0.9, mu = mu, size = size)
is_0.2 <- q0.9 - q0.1 + 5*(support < q0.1)*(q0.1 - support) + 5*(support > q0.9)*(support - q0.9)


# plot:
pdf("fig_intuition.pdf", width = 9, height = 2.8)

par(mfrow = c(1, 3), las = 1, mar = c(5, 5, 1, 4))

# logS:
plot(support, dens_F, xlab = "y", ylab = "", type = "l", ylim = c(0, 0.02), col = "white")
polygon(support, dens_F, col = col_F, border = col_F)
lsvals <- (-4):(-7)
axis(4, at = exp(lsvals), labels = lsvals)

par(las = 0)
mtext(expression(pred.~probability~p[y]), 2, 4, cex = 0.7)
mtext(expression(logS == log(p[y])), 4, 3, cex = 0.7)
abline(v = y, lty = "dashed")
arrows(x0 = y, x1 = max(support), y0 = dens_F[y + 1], y1 = dens_F[y + 1],
       length = 0.1, col = "steelblue")
text(x = y - 8, y = 0.0175, labels = "observed y", srt = 90)


# CRPS:
plot(support, cdf_F, type = "l", col = col_F2, lwd = 2, xlab = "y", ylab = "cumulative distr. function F(y)")

inds_above <- which(support >= y)
polygon(c(support[inds_above], y), c((1 - (1 - cdf_F[inds_above])^2), 1), col = "steelblue", border = NA)

inds_below <- which(support <= y)
polygon(c(support[inds_below], y), c(cdf_F[inds_below]^2, 0), col = "steelblue", border = NA)

lines(support, cdf_F^2, col = col_F2, lty = 4, lwd = 2)
lines(support, 1 - (1 - cdf_F)^2, col = col_F2, lty = 3, lwd = 2)

abline(v = y, lty = "dashed")
text(x = y - 8, y = 0.88, labels = "observed y", srt = 90)

legend("bottomright", legend = c(expression(1 - (1 - F(y))^2), "F(y)", expression(F(y)^2), "CRPS"),
       col = c(col_F2, col_F2, col_F2, "steelblue"), lty = c(3, 1, 4, NA), lwd = c(2, 2, 2, NA),
       pch = c(NA, NA, NA, 15),
       bty = "n")

# legend("bottomright", legend = c("F(y)", expression(F[y] %+-% F(y)(1 - F(y))), "CRPS"),
#        col = c(col_F2, col_F2, "steelblue"), lty = c(1, 3, NA), lwd = c(2, 2, NA),
#        pch = c(NA, NA, 15),
#        bty = "n")


# interval score:
plot(support, is_0.2, col = "steelblue", type = "l", lwd = 2, ylim = c(0, 600),
     xlab = "y", ylab = expression(IS[alpha]))
abline(v = c(q0.1, q0.9), col = col_F2, lty = 1, lwd = 2)
# text(105, 50, "}", cex = 5)
brackets(x1 = q0.9 + 5, x2 = q0.9 + 5, y2 = 0, y1 = q0.9 - q0.1, h = 5)
text(q0.9 + 25, 40, "u - l")
lines(c(150, 165, 165), c(is_0.2[151], is_0.2[151], is_0.2[166]))
text(175, 290, expression(slope == frac(2, alpha)))

text(q0.1 - 10, 300, "lower end l", srt = 90)
text(q0.9 - 10, 300, "upper end u", srt = 90)

dev.off()
