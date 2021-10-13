

plot(0, type = "n", xlim = c(-3, 3), ylim = c(0, 1), bty = "l",
     xlab = "Standardized effect size", ylab = "Density")
curve(2*dcauchy(x = x, location = 0, scale = sqrt(2)/2),
      from = 0, to = 3, lwd = 2,
      add = TRUE)
segments(-3, 0, 0, 0, lwd = 2)
segments(0, 0, 0, 2*dcauchy(0, 0, sqrt(2)/2), lwd = 2)
