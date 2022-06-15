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
data <- subset(data, select = c("study", "study_ID", "Lab", "Semantics", "Modality", "mean_age_1", "t", "n_1", "x_1", "x_2", "SD_1", "SD_2", "r"))
colnames(data)[1:6] <- c("study", "article", "lab", "semantics", "modality", "age")
data$semantics <- ordered(data$semantics, levels = c("Meaningless", "Meaningful"))
data$modality  <- ordered(data$modality)

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

# impute weighted correlations
# corr_missing <- is.na(data$corr_final)
# data$corr_final[corr_missing] <-
#   weighted.mean(data$corr_final[!corr_missing], 1/data$n_1[!corr_missing])

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
data$cohen_d_from_test_stat <- with(data, (t / sqrt(n_1)) * sqrt(2*(1-corr_final)))
data$cohen_d_final <- ifelse(
  !is.na(data$cohen_d_from_summary_stats),
  data$cohen_d_from_summary_stats,
  data$cohen_d_from_test_stat)

data$hedges_g <- with(data, cohen_d_final * (1 - 3 / (4*n_1 - 5)))
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

# The original random effect meta-analysis (https://osf.io/yu9c3/) -----
ran.ef.raw <- metafor::rma.mv(hedges_g, V = (hedges_g_se)^2, random = ~ article|lab,data = data)
ran.ef.raw


ran.ef.mod <- metafor::rma.mv(hedges_g, V = (hedges_g_se)^2, mods = ~ scale(age) * modality + scale(age) * semantics, random = ~ article|lab,data = data)
ran.ef.mod


# Limit the analysis only to semantics == 'meaningful' -----
data <- subset(data, semantics == "Meaningful")

## in JAGS -----
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
# Multivariate Meta-Analysis Model (k = 50; method: REML)
#
# Variance Components:
#
#             estim    sqrt  nlvls  fixed             factor
# sigma^2.1  0.0293  0.1713      8     no                lab
# sigma^2.2  0.0990  0.3146     12     no        lab/article
# sigma^2.3  0.0505  0.2247     50     no  lab/article/study
#
# Test for Heterogeneity:
#   Q(df = 49) = 183.3296, p-val < .0001
#
# Model Results:
#
#   estimate      se    zval    pval   ci.lb   ci.ub
#     0.3222  0.1228  2.6243  0.0087  0.0816  0.5628  **


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
# mu_pop          0.036  0.322   0.603 0.322 0.143   NA 0.003     2.0  2455 0.538 1.004
# sigma_lab       0.013  0.210   0.458 0.224 0.124   NA 0.002     1.9  2872 0.453 1.003
# sigma_article   0.119  0.306   0.544 0.317 0.108   NA 0.001     1.4  5322 0.262 1.001
# sigma_study     0.096  0.227   0.371 0.230 0.070   NA 0.001     1.9  2878 0.478 1.002

# Comparison:  median of Bayes <-> classical estimate:
# estimated effect size: 0.322 <-> 0.322, 95% CI: [0.036, 0.603] <-> [0.082, 0.563]
# sd of lab offset:      0.210 <-> 0.171
# sd of article offset:  0.306 <-> 0.315
# sd of study offset:    0.227 <-> 0.225
# => We can reproduce the results of the classical analysis, including CIs for the meta-analytic estimate

# Meta-analysis of the correlations ----
df <- subset(data, select = c("study", "article", "lab", "corr_z", "corr_z_se"))
df <- na.omit(df)
metafor::rma.mv(corr_z, V = (corr_z_se)^2, random = ~ 1|lab/article/study, data = df)
# Multivariate Meta-Analysis Model (k = 50; method: REML)
#
# Variance Components:
#
#             estim    sqrt  nlvls  fixed             factor
# sigma^2.1  0.0018  0.0427      8     no                lab
# sigma^2.2  0.0113  0.1061     12     no        lab/article
# sigma^2.3  0.0717  0.2678     50     no  lab/article/study
#
# Test for Heterogeneity:
#   Q(df = 49) = 88.4400, p-val = 0.0005
#
# Model Results:
#
#   estimate      se     zval    pval   ci.lb   ci.ub
#     0.7148  0.0709  10.0852  <.0001  0.5759  0.8537  ***

# Convert the estimate on the transformed scale into a correlation scale
tanh(c(estimate = 0.7148, lower = 0.5759, upper = 0.8537))
# estimate    lower    upper
#    0.614    0.520    0.693

# this results in the following point estimate of the within-subjects effect size:
0.322 / sqrt(2*(1-0.614)) # 0.366

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
# mu_pop          0.494  0.708   0.914 0.709 0.105   NA 0.001     1.1  8307 0.147 1.001
# mu_pop_corr     0.474  0.610   0.735 0.606 0.067   NA 0.001     1.1  8164 0.153 1.001
# sigma_lab       0.008  0.147   0.349 0.163 0.098   NA 0.002     1.5  4256 0.344 1.002
# sigma_article   0.010  0.137   0.303 0.148 0.082   NA 0.001     1.6  3814 0.380 1.002
# sigma_study     0.112  0.255   0.401 0.256 0.073   NA 0.001     1.7  3660 0.408 1.004


# Estimating the within-subjects effect size ----

## Meta-analysis of the within subjects effect size ----
### Frequentist analysis ----
#### Original model specification
df <- subset(data, select = c("study", "article", "lab", "within_g", "within_g_se"))
df <- na.omit(df)
metafor::rma.mv(within_g, V = (within_g_se)^2, random = ~ article|lab,data = df)
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.5078  0.1480  3.4302  0.0006  0.2176  0.7979  ***

#### Nested model
metafor::rma.mv(within_g, V = (within_g_se)^2, random = ~ 1|lab/article/study, data = df)
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.4602  0.1849  2.4884  0.0128  0.0977  0.8226  *

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
# mu_pop          0.014  0.446   0.836 0.443 0.207   NA 0.006     3.0  1121 0.752 1.003
# sigma_lab       0.008  0.206   0.511 0.233 0.145   NA 0.003     2.1  2167 0.548 1.001
# sigma_article   0.308  0.541   0.841 0.556 0.140   NA 0.002     1.2  7433 0.188 1.001
# sigma_study     0.263  0.384   0.528 0.389 0.068   NA 0.001     1.1  7938 0.170 1.000

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
# 19.767  0.443  0.196

hist(within_d, breaks = 100, freq = FALSE)
curve(dscaledt(x = x, df = o$par[[1]], center = o$par[[2]], scale = o$par[[3]]), add = TRUE, lwd = 5)

