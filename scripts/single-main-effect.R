library(BFDA)

# Global variables ----
# The following example assumes an informed prior: Normal(0.3, 0.1)
set.seed(2021)
mu_d <- 0.3
sd_d <- 0.1
sim_d <- rnorm(1e6, mu_d, sd_d)

# Using BFDA package ----
## fixed N ----
bfda_fixed_h0 <- BFDA.sim(expected.ES     = 0,
                              type        = "t.paired",
                              prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                              design      = "fixed.n",
                              n.max       = 200,
                              verbose     = TRUE)
BFDA.analyze(bfda_fixed_h0, design = "fixed", n = 200, boundary = 10)

bfda_fixed_h1 <- BFDA.sim(expected.ES    = sim_d,
                             type        = "t.paired",
                             prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                             design      = "fixed.n",
                             n.max       = 200,
                             verbose     = TRUE)
BFDA.analyze(bfda_fixed_h1, design = "fixed", n = 200, boundary = 10)
BFDA::evDens(BFDA.H1 = bfda_fixed_h1, BFDA.H0 = bfda_fixed_h0, n = 200, boundary = 10)

## sequential ----
bfda_seq_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                        design      = "sequential",
                        boundary    = 10,
                        n.min       = 20,
                        n.max       = 200,
                        stepsize    = 1,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq_h0, design = "sequential", n.min = 20, n.max = 200, boundary = 10)

bfda_seq_h1 <- BFDA.sim(expected.ES = sim_d,
                        type        = "t.paired",
                        prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                        design      = "sequential",
                        boundary    = 10,
                        n.min       = 20,
                        n.max       = 200,
                        stepsize    = 1,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq_h1, design = "sequential", n.min = 20, n.max = 200, boundary = 10)
plotBFDA(bfda_seq_h1, boundary = 10, n.min = 20, n.max = 200)

## sequential stepsize 20 ----
bfda_seq20_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                        design      = "sequential",
                        boundary    = 10,
                        n.min       = 20,
                        n.max       = 200,
                        stepsize    = 20,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq20_h0, design = "sequential", n.min = 20, n.max = 200, boundary = 10)

bfda_seq20_h1 <- BFDA.sim(expected.ES = sim_d,
                          type        = "t.paired",
                          prior       = list("normal", list(prior.mean = mu_d, prior.variance = sd_d^2)),
                          design      = "sequential",
                          boundary    = 10,
                          n.min       = 20,
                          n.max       = 200,
                          stepsize    = 20,
                          verbose     = TRUE,
                          cores       = 4,
                          ETA         = TRUE)
BFDA.analyze(bfda_seq20_h1, design = "sequential", n.min = 20, n.max = 200, boundary = 10)
plotBFDA(bfda_seq20_h1, boundary = 10, n.min = 20, n.max = 200)

## sequential default ----
bfda_seqdef_h0 <- BFDA.sim(expected.ES = 0,
                           type        = "t.paired",
                           design      = "sequential",
                           boundary    = 10,
                           n.min       = 20,
                           n.max       = 200,
                           stepsize    = 1,
                           verbose     = TRUE,
                           cores       = 4,
                           ETA         = TRUE)
BFDA.analyze(bfda_seqdef_h0, design = "sequential", n.min = 20, n.max = 200, boundary = 10)

bfda_seqdef_h1 <- BFDA.sim(expected.ES = sim_d,
                           type        = "t.paired",
                           design      = "sequential",
                           boundary    = 10,
                           n.min       = 20,
                           n.max       = 200,
                           stepsize    = 1,
                           verbose     = TRUE,
                           cores       = 4,
                           ETA         = TRUE)
BFDA.analyze(bfda_seqdef_h1, design = "sequential", n.min = 20, n.max = 200, boundary = 10)
plotBFDA(bfda_seqdef_h1, boundary = 10, n.min = 20, n.max = 200)

# Using BRMS package----
