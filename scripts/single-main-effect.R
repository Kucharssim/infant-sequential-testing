library(BFDA)

# Global variables ----
# The following example assumes an uninformed "default" prior, but one-sided
set.seed(2021)
prior_location <- 0
prior_scale    <- sqrt(2)/2
prior_alternative <- list("Cauchy", list(prior.location = prior_location, prior.scale = prior_scale))

d_of_interest  <- 0.5

min_n <- 5
max_n <- 50
step_size <- 5

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
BFDA.analyze(bfda_fixed_h1, design = "fixed", n = max_n, boundary = evidence_boundary)
BFDA::evDens(BFDA.H1 = bfda_fixed_h1, BFDA.H0 = bfda_fixed_h0, n = max_n, boundary = evidence_boundary)

## fixed N calculate desired sample size ----
desired_percentage <- 0.8
current_percentage <- 0
current_n <- max_n # we already know that for max_n we have less "power", so we start there

while(current_percentage < desired_percentage) {
  current_n <- current_n + 1
  bfda <- BFDA.sim(expected.ES = d_of_interest,
                   type        = "t.paired",
                   prior       = prior_alternative,
                   alternative = "greater",
                   design      = "fixed.n",
                   n.max       = current_n,
                   verbose     = FALSE,
                   cores       = 4)

  current_percentage <- mean(exp(bfda$sim$logBF) > evidence_boundary)
  cat(sprintf("Percentage of studies that showed evidence for H1 (BF > %s) for sample size %s: %s%%.\n",
      evidence_boundary, current_n, current_percentage*100)
      )
}

current_n
BFDA.analyze(bfda, design = "fixed", n = current_n, boundary = evidence_boundary)


bfda_0 <- BFDA.sim(expected.ES = 0,
                   type        = "t.paired",
                   prior       = prior_alternative,
                   alternative = "greater",
                   design      = "fixed.n",
                   n.max       = current_n)
BFDA.analyze(bfda_0, design = "fixed", n = current_n, boundary = evidence_boundary)

## sequential ----
bfda_seq_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        alternative = "greater",
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        cores       = 4)

bfda_seq_h1 <- BFDA.sim(expected.ES = d_of_interest,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        alternative = "greater",
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        cores       = 4)
BFDA.analyze(bfda_seq_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
BFDA.analyze(bfda_seq_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
plotBFDA(bfda_seq_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)

## sequential with stepsize ----
bfda_seqstep_h0 <- BFDA.sim(expected.ES = 0,
                            type        = "t.paired",
                            prior       = prior_alternative,
                            alternative = "greater",
                            design      = "sequential",
                            boundary    = evidence_boundary,
                            n.min       = min_n,
                            n.max       = max_n,
                            stepsize    = step_size,
                            cores       = 4)

bfda_seqstep_h1 <- BFDA.sim(expected.ES = d_of_interest,
                            type        = "t.paired",
                            prior       = prior_alternative,
                            alternative = "greater",
                            design      = "sequential",
                            boundary    = evidence_boundary,
                            n.min       = min_n,
                            n.max       = max_n,
                            stepsize    = step_size,
                            cores       = 4)
BFDA.analyze(bfda_seqstep_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
BFDA.analyze(bfda_seqstep_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
plotBFDA(bfda_seqstep_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)


