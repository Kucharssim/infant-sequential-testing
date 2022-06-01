library(readxl)
library(here)
library(runjags)
library(readxl)
library(metafor)

# Copy: https://osf.io/keq5a/ to your personal G Drive and export it as .xlsx file
data <- readxl::read_xlsx(
  path  = here::here("data", "RuleLearning MA Template.xlsx"),
  sheet = 1,
  na    = "NA")

# Convert F to t with an appropriate sign
data$t <- ifelse(
  !is.na(data$t),
  data$t,
  sqrt(abs(data$F)) * sign(data$x_1 - data$x_2)
)

data$study <- seq_len(nrow(data))
data <- subset(data, select = c("study", "study_ID", "Lab", "t", "n_1", "x_1", "x_2", "SD_1", "SD_2", "r"))
colnames(data)[1:3] <- c("study", "article", "lab")

# Calculate correlation using descriptive statistics and the within subject t-statistic
# Csibra, Hernik, Mascaro, Tatone, and Lengyel (2016)
calculate_corr <- function(n, t, m1, m2, sd1, sd2) {
  sum_of_variances <- sd1^2 + sd2^2
  inner_fraction <- n*(m1-m2)^2/(t^2)
  denominator <- 2*sd1*sd2

  result <- (sum_of_variances - inner_fraction) / denominator

  return(result)
}
data$corr_calculated <- with(data, calculate_corr(n_1, t, x_1, x_2, SD_1, SD_2))
plot(data$r, data$corr_calculated, pch = 21, cex = 1, bg = "grey")
abline(0, 1)

# if available, we take the reported correlation. Otherwise, we take the calculated correlation
data$corr_final     <- ifelse(!is.na(data$r), data$r, data$corr_calculated)

# Fisher's Z transformation of the correlation coefficients
data$corr_z <- atanh(data$corr_final)
data$corr_z_se <- 1 / sqrt(data$n_1 - 3)

# Calculate between subjects effect sizes
cohen_d <- function(m1, m2, sd1, sd2) {
  (m1 - m2)/pool_sd(sd1, sd2)
}
pool_sd <- function(sd1, sd2) {
  sqrt((sd1^2 + sd2^2)/2)
}

data$cohen_d_from_summary_stats <- with(data, cohen_d(x_1, x_2, SD_1, SD_2))
data$hedges_g <- with(data, cohen_d_from_summary_stats * (1 - 3 / (4*n_1 - 5)))
data$hedges_g_se <- with(data, sqrt(2*(1-corr_final)/n_1 + hedges_g^2/(2*n_1)))


# Calculate within subjects effect sizes
data$within_d_from_t_stat <- with(data, t / sqrt(n_1))
data$within_g <- with(
  data,
  ifelse(
    is.na(hedges_g),
    within_d_from_t_stat * (1 - 3 / (4*n_1 - 5)),
    hedges_g / sqrt(2*(1-corr_final))
  )
)
data$within_g_se <- with(
  data,
  sqrt((1/n_1 + within_g^2 / (2*n_1)) * 2*(1 - corr_final)) # Borenstein, Hedges, Higgins, & Rothstein, eq. 4.28
)

# Reproduce the original random effect meta-analysis (https://osf.io/yu9c3/) -----
ran.ef.raw <- metafor::rma.mv(hedges_g, V = (hedges_g_se)^2, random = ~ article|lab,data = data)
ran.ef.raw

## Reproduce in JAGS -----
# The correlation structure could be improved since the original analysis
# Assumed that the underlying effect size for a specific study is the same
# (using arguments random = ~ study|lab, struct = "CS")
# An alternative is a nested structure, i.e.,
# observed_effect_size[i] ~ dnorm(mu + offset_lab[l] + offset_article[a] + offset_study[i], se[i])
# offset_lab[l]     ~ dnorm(0, sigma_lab)
# offset_article[a] ~ dnorm(0, sigma_article)
# offset_study[i]   ~ dnorm(0, sigma_study)
df <- subset(data, select = c("study", "article", "lab", "hedges_g", "hedges_g_se"))
df <- na.omit(df)
metafor::rma.mv(hedges_g, V = (hedges_g_se)^2, random = ~ 1|lab/article/study, data = df)
# Variance Components:
#
#             estim    sqrt  nlvls  fixed             factor
# sigma^2.1  0.0286  0.1690      8     no                lab
# sigma^2.2  0.0475  0.2180     14     no        lab/article
# sigma^2.3  0.1054  0.3246     78     no  lab/article/study
#
# Model Results:
#
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.2504  0.0999  2.5069  0.0122  0.0546  0.4461  *


model <-
"
model{
  for(i in 1:n_lab) {
    lab_offset[i] ~ dnorm(0, tau_lab)
  }

  for(i in 1:n_article) {
    article_offset[i] ~ dnorm(0, tau_article)
  }

  for(i in 1:n_all) {
    study_offset[i] ~ dnorm(0, tau_study)
    d[i] ~ dnorm(mu_pop + lab_offset[lab[i]] + article_offset[article[i]] + study_offset[i], pow(se[i], -2))
  }

  mu_pop        ~ dnorm(0, 1)
  # using weakly informative priors for the variance components:
  # 1) It is unlikely that the variance is >> 1,
  #    more likely somewhere in the 0-1 range
  #    since the scale is a standardized effect size metric
  # 2) Since there is a relatively small number of labs and articles,
  #    the parameters are only weakly constrained by the data.
  #    Using 'flat' priors would drag the estimates to unreasonably high values.
  sigma_lab     ~ dgamma(2, 6)
  sigma_article ~ dgamma(2, 6)
  sigma_study   ~ dgamma(2, 6)

  tau_lab     <- pow(sigma_lab, -2)
  tau_article <- pow(sigma_article, -2)
  tau_study   <- pow(sigma_study, -2)
}
"

jags_data <- list(
  n_lab = length(unique(df$lab)),
  lab = as.integer(as.factor(df$lab)),
  n_article = length(unique(df$article)),
  article = as.integer(as.factor(df$article)),
  n_all = nrow(df),
  d = df$hedges_g,
  se = df$hedges_g_se
)

samples_between_subjects <- runjags::run.jags(
  model    = model,
  monitor  = c("mu_pop", "sigma_lab", "sigma_article", "sigma_study"),
  data     = jags_data,
  n.chains = 8,
  sample   = 10000
  )
summary(samples_between_subjects)
#               Lower95 Median Upper95  Mean    SD Mode MCerr MC%ofSD SSeff AC.10  psrf
# mu_pop          0.025  0.251   0.492 0.252 0.118   NA 0.002     2.1  2320 0.560 1.002
# sigma_lab       0.014  0.194   0.428 0.209 0.115   NA 0.002     1.8  3066 0.436 1.001
# sigma_article   0.043  0.224   0.431 0.231 0.101   NA 0.002     1.8  2978 0.449 1.004
# sigma_study     0.232  0.325   0.430 0.328 0.051   NA 0.001     1.1  7856 0.154 1.001

# Comparison:  median of Bayes <-> classical estimate:
# estimated effect size: 0.251 <-> 0.250, 95% CI: [0.025, 0.492] <-> [0.055, 0.446]
# sd of lab offset:      0.194 <-> 0.169
# sd of article offset:  0.224 <-> 0.218
# sd of study offset:    0.325 <-> 0.325
# => We can reproduce the results of the classical analysis, including CIs for the meta-analytic estimate

# Meta-analysis of the correlations ----
df <- subset(data, select = c("study", "article", "lab", "corr_z", "corr_z_se"))
df <- na.omit(df)
metafor::rma.mv(corr_z, V = (corr_z_se)^2, random = ~ 1|lab/article/study, data = df)
# Multivariate Meta-Analysis Model (k = 78; method: REML)
#
# Variance Components:
#
#             estim    sqrt  nlvls  fixed             factor
# sigma^2.1  0.0000  0.0000      8     no                lab
# sigma^2.2  0.0224  0.1497     14     no        lab/article
# sigma^2.3  0.0610  0.2470     78     no  lab/article/study
#
# Test for Heterogeneity:
#   Q(df = 77) = 132.6103, p-val < .0001
#
# Model Results:
#
# estimate      se     zval    pval   ci.lb   ci.ub     ​
#   0.7247  0.0636  11.4035  <.0001  0.6001  0.8493  ***

# Convert the estimate on the transformed scale into a correlation scale
tanh(c(estimate = 0.7247, lower = 0.6001, upper = 0.8493))
# estimate    lower    upper
#    0.620    0.537    0.691

# this results in the following point estimate of the within-subjects effect size:
0.25 / sqrt(2*(1-0.62)) # 0.29

model <-
"
model{
  for(i in 1:n_lab) {
    lab_offset[i] ~ dnorm(0, tau_lab)
  }

  for(i in 1:n_article) {
    article_offset[i] ~ dnorm(0, tau_article)
  }

  for(i in 1:n_all) {
    study_offset[i] ~ dnorm(0, tau_study)
    z[i] ~ dnorm(mu_pop + lab_offset[lab[i]] + article_offset[article[i]] + study_offset[i], pow(se[i], -2))
  }

  mu_pop        ~ dnorm(0, 1)
  mu_pop_corr   <- tanh(mu_pop)

  sigma_lab     ~ dgamma(2, 6)
  sigma_article ~ dgamma(2, 6)
  sigma_study   ~ dgamma(2, 6)

  tau_lab     <- pow(sigma_lab, -2)
  tau_article <- pow(sigma_article, -2)
  tau_study   <- pow(sigma_study, -2)
}
"

df <- data
jags_data <- list(
  n_lab     = length(unique(df$lab)),
  lab       = as.integer(as.factor(df$lab)),
  n_article = length(unique(df$article)),
  article   = as.integer(as.factor(df$article)),
  n_all     = nrow(df),
  z         = df$corr_z,
  se        = df$corr_z_se
)

samples_correlation <- runjags::run.jags(
  model    = model,
  monitor  = c("mu_pop", "mu_pop_corr",  "sigma_lab", "sigma_article", "sigma_study"),
  data     = jags_data,
  n.chains = 8,
  sample   = 10000
)
summary(samples_correlation)
#               Lower95 Median Upper95  Mean    SD Mode MCerr MC%ofSD SSeff AC.10  psrf
# mu_pop          0.551  0.725   0.908 0.727 0.089   NA 0.001     1.2  7378 0.167 1.001
# mu_pop_corr     0.511  0.620   0.729 0.618 0.055   NA 0.001     1.2  7511 0.159 1.001
# sigma_lab       0.006  0.123   0.297 0.138 0.083   NA 0.001     1.5  4525 0.331 1.001
# sigma_article   0.021  0.154   0.302 0.160 0.075   NA 0.001     1.6  4116 0.361 1.001
# sigma_study     0.129  0.244   0.356 0.244 0.058   NA 0.001     1.7  3411 0.439 1.001


# Estimating the within-subjects effect size ----

## Meta-analysis of the within subjects effect size ----
### Frequentist analysis ----
#### Original model specification
df <- subset(data, select = c("study", "article", "lab", "within_g", "within_g_se"))
df <- na.omit(df)
metafor::rma.mv(within_g, V = (within_g_se)^2, random = ~ article|lab,data = df)
# estimate      se    zval    pval   ci.lb   ci.ub    ​
#   0.4197  0.1313  3.1960  0.0014  0.1623  0.6770  **

#### Nested model
metafor::rma.mv(within_g, V = (within_g_se)^2, random = ~ 1|lab/article/study, data = df)
# estimate      se    zval    pval   ci.lb   ci.ub   ​
#   0.3574  0.1389  2.5727  0.0101  0.0851  0.6297  *

### Bayesian analysis ----

model <-
"
model{
  for(i in 1:n_lab) {
    lab_offset[i] ~ dnorm(0, tau_lab)
  }

  for(i in 1:n_article) {
    article_offset[i] ~ dnorm(0, tau_article)
  }

  for(i in 1:n_all) {
    study_offset[i] ~ dnorm(0, tau_study)
    d[i] ~ dnorm(mu_pop + lab_offset[lab[i]] + article_offset[article[i]] + study_offset[i], pow(se[i], -2))
  }

  mu_pop        ~ dnorm(0, 1)
  # using weakly informative priors for the variance components:
  # 1) It is unlikely that the variance is >> 1,
  #    more likely somewhere in the 0-1 range
  #    since the scale is a standardized effect size metric
  # 2) Since there is a relatively small number of labs and articles,
  #    the parameters are only weakly constrained by the data.
  #    Using 'flat' priors would drag the estimates to unreasonably high values.
  sigma_lab     ~ dgamma(2, 6)
  sigma_article ~ dgamma(2, 6)
  sigma_study   ~ dgamma(2, 6)

  tau_lab     <- pow(sigma_lab, -2)
  tau_article <- pow(sigma_article, -2)
  tau_study   <- pow(sigma_study, -2)
}
"

jags_data <- list(
  n_lab = length(unique(df$lab)),
  lab = as.integer(as.factor(df$lab)),
  n_article = length(unique(df$article)),
  article = as.integer(as.factor(df$article)),
  n_all = nrow(df),
  d = df$within_g,
  se = df$within_g_se
)

samples_within <- runjags::run.jags(
  model    = model,
  monitor  = c("mu_pop", "sigma_lab", "sigma_article", "sigma_study"),
  data     = jags_data,
  n.chains = 8,
  sample   = 10000
)
summary(samples_within)
#               Lower95 Median Upper95  Mean    SD Mode MCerr MC%ofSD SSeff AC.10  psrf
# mu_pop          0.031  0.359   0.688 0.359 0.165   NA 0.005     2.9  1185 0.746 1.008
# sigma_lab       0.010  0.200   0.487 0.224 0.138   NA 0.003     2.2  2050 0.553 1.006
# sigma_article   0.143  0.390   0.673 0.399 0.133   NA 0.002     1.9  2896 0.429 1.003
# sigma_study     0.400  0.509   0.631 0.513 0.059   NA 0.001     1.1  8486 0.140 1.001

within_d <- do.call(rbind.data.frame, samples_within$mcmc)$mu_pop
hist(within_d, breaks = 100, freq = FALSE)
lines(density(within_d), lwd = 5)

##### Finding an appropriate representation of the prior distribution ----
# Now we have an estimate of the effect size, but we need to represent it in terms
# of a probability distribution so that it can be used as a prior.
# Here, we fit a scaled-shifted t-distribution to the MCMC samples of the
# within-subjects effect size.
dscaledt <- function(x, df, center, scale, log = FALSE) {
  out <- dt(x = (x-center)/scale, df = df, ncp = 0, log = FALSE) / scale

  if(log) {
    return(log(out))
  } else {
    return(out)
  }
}

fn <- function(pars, data) {
  df     <- pars[[1]]
  center <- pars[[2]]
  scale  <- pars[[3]]

  ll <- dscaledt(x = data, df = df, center = center, scale = scale, log = TRUE)

  return(sum(ll))
}

o <- optim(par = list(df = 2, center = 0, scale = 1), fn = fn, data = within_d,
           lower = c(1, -Inf, 0.001), upper = c(Inf, Inf, Inf),
           method = "L-BFGS-B", control = list(fnscale = -1))
#     df center  scale
# 10.293  0.358  0.148

hist(within_d, breaks = 100, freq = FALSE)
curve(dscaledt(x = x, df = o$par[[1]], center = o$par[[2]], scale = o$par[[3]]), add = TRUE, lwd = 5)

