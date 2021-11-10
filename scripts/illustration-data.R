library(BayesFactor)
source("R/helpers.R") # loads functions getBF() and sequential_plot()

data <- read.csv("data/data_example.csv")
head(data)


# Fixed-n design
fixed <- getBF(data)

# Sequential design
sequential <- data.frame(
  sampleSize = 1:nrow(data),
  bf         = rep(1, nrow(data))
)
for(i in 2:nrow(data)) {
  sequential$bf[i] <- getBF(data[seq_len(i),])
}

sequential_plot()
points(sequential$sampleSize, log(sequential$bf), pch = 21, bg = ifelse(sequential$bf > 1/10, "gray", "red"), cex = 1, lwd = 1.3)


# Sequential design in batches of five
sequential_batches <- data.frame(
  sampleSize = seq(5, nrow(data), by = 5),
  bf         = 1
)
for(i in 1:nrow(sequential_batches)) {
  sequential_batches$bf[i] <- getBF(data[seq_len(sequential_batches$sampleSize[i]),])
}

sequential_plot()
points(sequential_batches$sampleSize, log(sequential_batches$bf), pch = 21, bg = ifelse(sequential_batches$bf > 1/10, "gray", "red"), cex = 1, lwd = 1.3)
