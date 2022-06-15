library(BFDA)
sink(file = "default-prior.Rout")

# Global variables ----
# The following example assumes an uninformed "default" prior, but one-sided
set.seed(2021)
prior_location <- 0
prior_scale    <- sqrt(2)/2
prior_alternative <- list("Cauchy", list(prior.location = prior_location, prior.scale = prior_scale))

d_of_interest  <- 0.45

min_n <- 5
max_n <- 50
step_size <- 5

evidence_boundary <- 10
verbose <- FALSE
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
                   verbose     = verbose,
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

# Add nice plots ----
library(ggplot2)
library(dplyr)

## Figure 3: Fixed N, distribution of Bayes factors ----
df <- rbind(
  bfda_fixed_h0$sim,
  bfda_fixed_h1$sim
)
df$evidence <- case_when(
  df$logBF < log(1/10) ~ "Null",
  df$logBF > log(10)   ~ "Alternative",
  TRUE ~ "Undecisive"
) |>
  as.factor()
df$trueHypothesis <- ifelse(df$true.ES == 0, "Under null", "Under alternative") |>
  factor(levels = c("Under null", "Under alternative"))


ggplot(df, aes(x = logBF, fill = evidence)) +
  geom_histogram(col = "black", boundary = log(1/10), binwidth = (log(10) - log(1/10))/10) +
  geom_rug() +
  geom_vline(xintercept = 0, linetype = 2, size = 0.8) +
  geom_vline(xintercept = log(c(1/10, 10)), linetype = 3, size = 0.8) +
  scale_x_continuous(
    limits = log(c(1e-2, 1e5)),
    breaks = log(c(1/100, 1/10, 1, 10, 1e2, 1e3, 1e4, 1e5)),
    labels = c(expression(10^-2), expression(10^-1), "1", expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5))
  ) +
  scale_fill_manual(
    name = "Evidence",
    values = c("steelblue", "gray", "yellow2"), # somewhat color-blind friendly setup
    breaks = c("Null", "Undecisive", "Alternative")
  ) +
  facet_grid(trueHypothesis~.) +
  theme_classic(base_size = 18) +
  ylab("Count") +
  xlab("log(BF)")
ggsave(filename = "single-main-effect_fixed-n_bf-histogram.png",
       path = here::here("figures"), width = 7, height = 5)


## Figure 4: Sequential design.
library(patchwork)

### Under null ----
df <- bfda_seq_h0$sim

df$logBF <- case_when(
  df$logBF < log(1/evidence_boundary) ~ log(1/evidence_boundary),
  df$logBF > log(evidence_boundary) ~ log(evidence_boundary),
  TRUE ~ df$logBF

)
df$crossed0 <- df$logBF <= log(1/10)
df$crossed1 <- df$logBF >= log(10)

xmax <- 60

accumulation <- ggplot(df, aes(x = n, y = logBF, group = id)) +
  geom_hline(yintercept = log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)), linetype = 3, size = 0.3) +
  geom_hline(yintercept = log(c(1/10, 10)), linetype = 2, size = 0.5) +
  geom_point(aes(x = n, y = log(1/10)), data = subset(df, crossed0)) +
  geom_point(aes(x = n, y = log(10)),   data = subset(df, crossed1)) +
  geom_point(aes(x = 50, y = logBF), data = subset(df, !crossed0 & !crossed1 & n == 50)) +
  geom_line(alpha = 0.05) +
  geom_text(
    aes(x = x, y = y, label = text),
    data = data.frame(
      x = xmax,
      y = c(
        mean(log(c(1/30, 1/10))),
        mean(log(c(1/10, 1/3))),
        mean(log(c(1/3,  1))),
        mean(log(c(1,    3))),
        mean(log(c(3,    10))),
        mean(log(c(10,   30)))
      ),
      text = c(
        "Strong H0",
        "Moderate H0",
        "Anecdotal H0",
        "Anecdotal H1",
        "Moderate H1",
        "Strong H1"
      )
    ),
    inherit.aes = FALSE,
    hjust = 0.8,
    size = 5.5
  ) +
  scale_x_continuous(
    name = NULL,
    guide = guide_none(),
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
  ) +
  scale_y_continuous(
    name = "log(BF)",
    breaks = log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)),
    labels = c("1/30", "1/10", "1/3", "1", "3", "10", "30")
  ) +
  theme_classic(base_size = 20)

accumulation

sampleSize0 <- ggplot(data = subset(df, crossed0), aes(x = n)) +
  geom_histogram(breaks = seq(0, 50, by = 2.5), col = "black", fill = "gray", closed = "left") +
  scale_x_continuous(
    name = "Sample size",
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
  ) +
  scale_y_reverse(
    name = NULL, guide = guide_none()
  ) +
  theme_classic(base_size = 20)
sampleSize0

sampleSize1 <- ggplot(data = subset(df, crossed1), aes(x = n)) +
  geom_histogram(breaks = seq(0, 50, by = 2.5), col = "black", fill = "gray", closed = "left") +
  scale_x_continuous(
    name = NULL,
    guide = guide_none(),
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
    ) +
  scale_y_continuous(
    name = NULL,
    guide = guide_none()
  ) +
  ggtitle("Under null") +
  theme_classic(base_size = 20)
sampleSize1

evidenceAtStop <- ggplot(data = subset(df, n == 50), aes(x = logBF)) +
  geom_histogram(
    breaks = seq(log(1/evidence_boundary), log(evidence_boundary), length.out = 10),
    col = "black", fill = "gray"
    ) +
  theme_classic(base_size = 18) +
  coord_flip() +
  scale_y_continuous(
    name = NULL,
    guide = guide_none()
  ) +
  scale_x_continuous(
    name   = NULL,
    guide  = guide_none(),
    limits = c(log(1/30), log(30))
  )

evidenceAtStop

sampleSize1 + plot_spacer() +
  accumulation +  evidenceAtStop +
  sampleSize0 + plot_spacer() +
  patchwork::plot_layout(ncol = 2, heights = c(0.15, 0.7, 0.15), widths = c(0.85, 0.15))# +
  #patchwork::plot_annotation(title = "Under null")

ggsave(filename = "single-main-effect_sequential-h0_reworked.png",
       path = here::here("figures"),
       width = 10, height = 7)


### Under alternative ----
df <- bfda_seq_h1$sim

df$logBF <- case_when(
  df$logBF < log(1/evidence_boundary) ~ log(1/evidence_boundary),
  df$logBF > log(evidence_boundary) ~ log(evidence_boundary),
  TRUE ~ df$logBF

)
df$crossed0 <- df$logBF <= log(1/10)
df$crossed1 <- df$logBF >= log(10)

xmax <- 60

accumulation <- ggplot(df, aes(x = n, y = logBF, group = id)) +
  geom_hline(yintercept = log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)), linetype = 3, size = 0.3) +
  geom_hline(yintercept = log(c(1/10, 10)), linetype = 2, size = 0.5) +
  #geom_point(aes(x = n, y = log(1/10)), data = subset(df, crossed0)) +
  geom_point(aes(x = n, y = log(10)),   data = subset(df, crossed1)) +
  geom_point(aes(x = 50, y = logBF), data = subset(df, !crossed0 & !crossed1 & n == 50)) +
  geom_line(alpha = 0.05) +
  geom_text(
    aes(x = x, y = y, label = text),
    data = data.frame(
      x = xmax,
      y = c(
        mean(log(c(1/30, 1/10))),
        mean(log(c(1/10, 1/3))),
        mean(log(c(1/3,  1))),
        mean(log(c(1,    3))),
        mean(log(c(3,    10))),
        mean(log(c(10,   30)))
      ),
      text = c(
        "Strong H0",
        "Moderate H0",
        "Anecdotal H0",
        "Anecdotal H1",
        "Moderate H1",
        "Strong H1"
      )
    ),
    inherit.aes = FALSE,
    hjust = 0.8,
    size = 5.5
  ) +
  scale_x_continuous(
    name = NULL,
    guide = guide_none(),
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
  ) +
  scale_y_continuous(
    name = "log(BF)",
    breaks = log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)),
    labels = c("1/30", "1/10", "1/3", "1", "3", "10", "30")
  ) +
  theme_classic(base_size = 20)

accumulation

sampleSize0 <- ggplot(data = subset(df, crossed0), aes(x = n)) +
  geom_histogram(breaks = seq(0, 50, by = 2.5), col = "black", fill = "gray", closed = "left") +
  scale_x_continuous(
    name = "Sample size",
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
  ) +
  scale_y_reverse(
    name = NULL, guide = guide_none(), limits = c(0, 0)
  ) +
  theme_classic(base_size = 20)
sampleSize0

sampleSize1 <- ggplot(data = subset(df, crossed1), aes(x = n)) +
  geom_histogram(breaks = seq(0, 50, by = 2.5), col = "black", fill = "gray", closed = "left") +
  scale_x_continuous(
    name = NULL,
    guide = guide_none(),
    breaks = seq(0, 50, by = 10),
    limits = c(0, xmax)
  ) +
  scale_y_continuous(
    name = NULL,
    guide = guide_none()
  ) +
  ggtitle("Under alternative") +
  theme_classic(base_size = 20)
sampleSize1

evidenceAtStop <- ggplot(data = subset(df, n == 50), aes(x = logBF)) +
  geom_histogram(
    breaks = seq(log(1/evidence_boundary), log(evidence_boundary), length.out = 10),
    col = "black", fill = "gray"
  ) +
  theme_classic(base_size = 18) +
  coord_flip() +
  scale_y_continuous(
    name = NULL,
    guide = guide_none()
  ) +
  scale_x_continuous(
    name   = NULL,
    guide  = guide_none(),
    limits = c(log(1/30), log(30))
  )

evidenceAtStop

sampleSize1 + plot_spacer() +
  accumulation +  evidenceAtStop +
  sampleSize0 + plot_spacer() +
  patchwork::plot_layout(ncol = 2, heights = c(0.15, 0.7, 0), widths = c(0.85, 0.15))# +
#patchwork::plot_annotation(title = "Under null")

ggsave(filename = "single-main-effect_sequential-h1_reworked.png",
       path = here::here("figures"),
       width = 10, height = 7)
