library(BFDA)

# Global variables ----
# The following example assumes an uninformed "default" prior, but one-sided
set.seed(2021)
prior_location <- 0
prior_scale    <- sqrt(2)/2
prior_alternative <- list("Cauchy", list(prior.location = prior_location, prior.scale = prior_scale))

d_of_interest  <- 0.5

min_n <- 20
max_n <- 100

evidence_boundary <- 10
# Using BFDA package ----
## fixed N ----
bfda_fixed_h0 <- BFDA.sim(expected.ES = 0,
                          type        = "t.paired",
                          prior       = prior_alternative,
                          alternative = "greater",
                          design      = "fixed.n",
                          n.max       = max_n,
                          verbose     = TRUE)

bfda_fixed_h1 <- BFDA.sim(expected.ES = d_of_interest,
                          type        = "t.paired",
                          prior       = prior_alternative,
                          alternative = "greater",
                          design      = "fixed.n",
                          n.max       = max_n,
                          verbose     = TRUE)
BFDA.analyze(bfda_fixed_h0, design = "fixed", n = max_n, boundary = evidence_boundary)
BFDA.analyze(bfda_fixed_h1, design = "fixed", n = nax_n, boundary = evidence_boundary)
BFDA::evDens(BFDA.H1 = bfda_fixed_h1, BFDA.H0 = bfda_fixed_h0, n = max_n, boundary = evidence_boundary)

## sequential ----
bfda_seq_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)

bfda_seq_h1 <- BFDA.sim(expected.ES = sim_d,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
plotBFDA(bfda_seq_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)

## sequential stepsize 20 ----
bfda_seq20_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 20,
                        verbose     = TRUE,
                        cores       = 4,
                        ETA         = TRUE)
BFDA.analyze(bfda_seq20_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)

bfda_seq20_h1 <- BFDA.sim(expected.ES = sim_d,
                          type        = "t.paired",
                          prior       = prior_alternative,
                          design      = "sequential",
                          boundary    = evidence_boundary,
                          n.min       = min_n,
                          n.max       = max_n,
                          stepsize    = 20,
                          verbose     = TRUE,
                          cores       = 4,
                          ETA         = TRUE)
BFDA.analyze(bfda_seq20_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
plotBFDA(bfda_seq20_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)


# Using BRMS package----
