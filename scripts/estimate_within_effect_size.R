library(readxl)
library(here)
library(runjags)

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
data$hedges_g_se <- with(data, sqrt(2*(1-data$corr_final)/n_1 + hedges_g^2/(2*n_1)))


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
## Assuming independence between the between subject effect size and the correlation ----
between_d <- do.call(rbind.data.frame, samples_between_subjects$mcmc)$mu_pop
corr      <- do.call(rbind.data.frame, samples_correlation$mcmc)$mu_pop_corr
within_d  <- between_d / (2*(1-corr))
hist(within_d)
c(quantile(within_d, c(0.025, 0.25, 0.5, 0.75, 0.975)), mean = mean(within_d), sd = sd(within_d))


## Meta-analysis of the within subjects effect size ----
# Here we estimate the within-subject effect size again but taking into account
# the correlation between the between-subjects effect size and the correlation
# between the two paired measurements.

model <-
"
model{
  for(i in 1:3) { # lab, article, study
    for(j in 1:2) { # corr, hedges g
      sigma[i, j] ~ dgamma(2, 6)
    }
    half_rho[i] ~ dbeta(2, 2)
    rho[i] = half_rho[i] * 2 - 1

    Sigma[i, 1, 1] = sigma[i, 1]*sigma[i, 1]
    Sigma[i, 1, 2] = sigma[i, 1]*sigma[i, 2]*rho[i]
    Sigma[i, 2, 1] = sigma[i, 2]*sigma[i, 1]*rho[i]
    Sigma[i, 2, 2] = sigma[i, 2]*sigma[i, 2]
  }

  for(i in 1:n_lab) {
    lab_offset[i,1:2] ~ dmnorm.vcov(c(0, 0), Sigma[1,1:2,1:2])
  }

  for(i in 1:n_article) {
    article_offset[i,1:2] ~ dmnorm.vcov(c(0, 0), Sigma[2,1:2,1:2])
  }

  for(i in 1:n_all) {
    study_offset[i,1:2] ~ dmnorm.vcov(c(0, 0), Sigma[3,1:2,1:2])

    corr_z[i] ~ dnorm(corr_mu_pop + lab_offset[lab[i],1] + article_offset[article[i],1] + study_offset[i,1], pow(corr_se[i], -2))
    corr[i] <- tanh(corr_z[i])

    betw_d[i] ~ dnorm(betw_mu_pop + lab_offset[lab[i],2] + article_offset[article[i],2] + study_offset[i,2], pow(betw_se[i], -2))
  }

  # priors for correlations
  corr_mu_pop        ~ dnorm(0, 1)
  corr_mu_pop_corr   <- tanh(corr_mu_pop)

  # priors for between subjects
  betw_mu_pop        ~ dnorm(0, 1)

  # calculate within subjects
  with_mu_pop <- betw_mu_pop / (2 * (1-corr_mu_pop_corr))
}
"
df <- subset(data, select = c("study", "article", "lab", "hedges_g", "hedges_g_se", "corr_z", "corr_z_se"))
df <- subset(df, subset = !is.na(hedges_g_se))
jags_data <- list(
  n_all     = nrow(df),
  n_lab     = length(unique(df$lab)),
  n_article = length(unique(df$article)),
  corr_z    = df$corr_z,
  corr_se   = df$corr_z_se,
  betw_d    = df$hedges_g,
  betw_se   = df$hedges_g_se,
  lab       = as.integer(as.factor(df$lab)),
  article   = as.integer(as.factor(df$article))
)

samples_within <- runjags::run.jags(
  model    = model,
  monitor  = c("corr_mu_pop_corr", "betw_mu_pop", "with_mu_pop", "rho", "sigma"),
  data     = jags_data,
  n.chains = 8,
  sample   = 20000,
  thin     = 2
)
summary(samples_within)
#                  Lower95 Median Upper95   Mean    SD Mode MCerr MC%ofSD SSeff AC.10  psrf
# corr_mu_pop_corr   0.507  0.620   0.726  0.617 0.056   NA 0.001     2.6  1538 0.618 1.005
# betw_mu_pop       -0.002  0.252   0.519  0.257 0.129   NA 0.005     4.0   624 0.832 1.019
# with_mu_pop       -0.008  0.328   0.723  0.344 0.189   NA 0.007     3.9   675 0.816 1.020
# rho[1]            -0.863 -0.107   0.708 -0.088 0.432   NA 0.008     1.7  3269 0.364 1.001
# rho[2]            -0.656  0.164   0.859  0.137 0.409   NA 0.009     2.1  2223 0.503 1.002
# rho[3]            -0.201  0.245   0.682  0.242 0.226   NA 0.006     2.5  1566 0.604 1.003
# sigma[1,1]         0.006  0.124   0.297  0.138 0.084   NA 0.003     3.2   994 0.626 1.004
# sigma[2,1]         0.018  0.149   0.296  0.156 0.075   NA 0.002     3.2   957 0.670 1.005
# sigma[3,1]         0.120  0.243   0.354  0.242 0.059   NA 0.002     3.5   800 0.763 1.005
# sigma[1,2]         0.013  0.206   0.455  0.223 0.123   NA 0.005     4.1   588 0.651 1.036
# sigma[2,2]         0.028  0.218   0.422  0.227 0.103   NA 0.005     4.4   516 0.716 1.010
# sigma[3,2]         0.227  0.326   0.428  0.328 0.051   NA 0.001     2.3  1902 0.530 1.003
within_d <- do.call(rbind.data.frame, samples_within$mcmc)$with_mu_pop
hist(within_d, breaks = 100, freq = FALSE)


### Finding an appropriate representation of the prior distribution ----
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
#    df center  scale
# 5.646  0.333  0.151

hist(within_d, breaks = 100, freq = FALSE)
curve(dscaledt(x = x, df = o$par[[1]], center = o$par[[2]], scale = o$par[[3]]), add = TRUE, lwd = 2)

