library(BFDA)
sink(file = "informative-prior.Rout")
# Global variables ----
# The following example assumes an uninformed "default" prior, but one-sided
set.seed(2021)
# the prior is determined in the script estimate_within_effect_size.R
prior_alternative <- list("t", list(prior.df = 19.767, prior.location = 0.443, prior.scale = 0.196))

d_of_interest  <- 0.45

min_n <- 5
max_n <- 50
step_size <- 5

evidence_boundary <- 10
verbose = FALSE
# Using BFDA package ----
## fixed N ----
cat("Fixed N design\n\n")
bfda_fixed_h0 <- BFDA.sim(expected.ES = 0,
                          type        = "t.paired",
                          prior       = prior_alternative,
                          alternative = "greater",
                          design      = "fixed.n",
                          n.max       = max_n,
                          verbose     = verbose)

bfda_fixed_h1 <- BFDA.sim(expected.ES = d_of_interest,
                          type        = "t.paired",
                          prior       = prior_alternative,
                          alternative = "greater",
                          design      = "fixed.n",
                          n.max       = max_n,
                          verbose     = verbose)
cat("Fixed N design under null \n")
BFDA.analyze(bfda_fixed_h0, design = "fixed", n = max_n, boundary = evidence_boundary)
cat("Fixed N design under alternative \n")
BFDA.analyze(bfda_fixed_h1, design = "fixed", n = max_n, boundary = evidence_boundary)
BFDA::evDens(BFDA.H1 = bfda_fixed_h1, BFDA.H0 = bfda_fixed_h0, n = max_n, boundary = evidence_boundary)

## fixed N calculate desired sample size ----
cat("Fixed N design, calculate desired sample size \n\n")
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
                   n.max       = current_n,
                   verbose     = verbose)
BFDA.analyze(bfda_0, design = "fixed", n = current_n, boundary = evidence_boundary)

## sequential ----
cat("Sequential design \n\n")
bfda_seq_h0 <- BFDA.sim(expected.ES = 0,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        alternative = "greater",
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        cores       = 4,
                        verbose     = verbose)

bfda_seq_h1 <- BFDA.sim(expected.ES = d_of_interest,
                        type        = "t.paired",
                        prior       = prior_alternative,
                        alternative = "greater",
                        design      = "sequential",
                        boundary    = evidence_boundary,
                        n.min       = min_n,
                        n.max       = max_n,
                        stepsize    = 1,
                        cores       = 4,
                        verbose     = verbose)
cat("Sequential design under null \n")
BFDA.analyze(bfda_seq_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
cat("Sequential design under alternative \n")
BFDA.analyze(bfda_seq_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)


plotBFDA(bfda_seq_h0, boundary = evidence_boundary, n.min = min_n, n.max = max_n)
mtext(expression("Under"~H[0]))
plotBFDA(bfda_seq_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)
mtext(expression("Under"~H[1]))

## sequential with stepsize ----
cat("Sequential design with steps \n\n")
bfda_seqstep_h0 <- BFDA.sim(expected.ES = 0,
                            type        = "t.paired",
                            prior       = prior_alternative,
                            alternative = "greater",
                            design      = "sequential",
                            boundary    = evidence_boundary,
                            n.min       = min_n,
                            n.max       = max_n,
                            stepsize    = step_size,
                            cores       = 4,
                            verbose     = verbose)

bfda_seqstep_h1 <- BFDA.sim(expected.ES = d_of_interest,
                            type        = "t.paired",
                            prior       = prior_alternative,
                            alternative = "greater",
                            design      = "sequential",
                            boundary    = evidence_boundary,
                            n.min       = min_n,
                            n.max       = max_n,
                            stepsize    = step_size,
                            cores       = 4,
                            verbose     = verbose)
cat("Sequential design with steps under null \n")
BFDA.analyze(bfda_seqstep_h0, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
cat("Sequential design with steps under alternative \n")
BFDA.analyze(bfda_seqstep_h1, design = "sequential", n.min = min_n, n.max = max_n, boundary = evidence_boundary)
plotBFDA(bfda_seqstep_h1, boundary = evidence_boundary, n.min = min_n, n.max = max_n)


sink()
