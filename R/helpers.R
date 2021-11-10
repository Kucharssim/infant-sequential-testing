getBF <- function(data) {
  # first we make a simple helper function that extracts the relevant BF
  result <- ttestBF(
    x = data$inconsistent,
    y = data$consistent,
    paired = TRUE,
    rscale = "medium", # this corresponds to the Cauchy(0, sqrt(2)/2) prior
    nullInterval = c(0, Inf) # setting one-sided test
  )
  bf <- extractBF(result, onlybf = TRUE)[1]

  return(bf)
}

sequential_plot <- function() {
  # adapted from https://www.shinyapps.org/apps/RGraphCompendium/index.php#evidential-flow
  # for this specific analysis
  par(cex = 1.8, mar = c(4, 4, 2, 6))
  plot(x = sequential$sampleSize,
       y = log(sequential$bf),
       ylim = log(c(1/30, 10)),
       bty = "n", axes = FALSE,
       ylab = "",
       xlab = "",
       type = "n")
  evidence_bounds <- log(c(1/30, 1/10, 1/3, 1, 3, 10))
  axis(side = 1, at = seq(0, 40, 10))
  axis(side = 2, at = evidence_bounds, labels = c("1/30", "1/10", "1/3", "1", "3", "10"), las = 1)
  axis(side = 4, at = evidence_bounds, labels = rep("", 6))
  axis(side = 4, at = sapply(2:length(evidence_bounds), function(i) mean(evidence_bounds[(i-1):i])), tick = FALSE,
       labels = c("Strong", "Moderate", "Anecdotal", "Anecdotal", "Strong"), las = 1, line = -0.8, cex = 1.5)

  mtext("n", side = 1, line = 2, cex = 2)
  mtext(expression(BF[10]), side = 2, line = 2.5, cex = 2)
  mtext("Evidence", side = 4, line = 4.5, cex = 2)

  for(evid in evidence_bounds) segments(0, evid, 40, evid, lty = 2, col = "gray")
  segments(0, 0, 40, 0, lwd = 2)
  segments(0, log(1/10), 40, log(1/10), lty = 2)
  segments(0, log(10),   40, log(10),   lty = 2)

  arrows(5, log(1/12), 5, log(1/26), length = 0.15, lwd = 2.5)
  text(12.5, mean(c(log(1/12), log(1/26))), expression("Evidence for"~H[0]), cex = 0.9)

  arrows(5, log(4), 5, log(8), length = 0.15, lwd = 2.5)
  text(12.5, mean(c(log(4), log(8))), expression("Evidence for"~H[1]), cex = 0.9)
}
