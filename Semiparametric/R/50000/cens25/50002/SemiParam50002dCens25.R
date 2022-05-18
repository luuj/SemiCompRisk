#####################################################################################
#####################################################################################
# CausalSemiComp
# Choose scenario in line 13. Scenarios are defined in the GetScenarioParams function
#####################################################################################
rm(list = ls())
library(dplyr)
library(survival)
library(CausalSemiComp)
library(boot)

lett <- "d"
scen <- 50002
set.seed(scen + match(lett,letters))

#####################################################################################
params <- GetScenarioParams(scenario.num = scen)

no.protected <- F
no.large <- F
cens.rate <- 0.026
n.sample <- 2000
max.iter <- 5000
n.sim <- 100
H.times <- seq(0.25, 25, 0.25) # values at which all 6 baseline hazard functions are estimated + bootstrapped
tau <- 15
R <- 100
n.gamma.vals <- 100
n.sample.pers <- 100
X.x1 <- rep(seq(-2, 2, 0.5), times = 2)
X.x2 <- rep(0:1, each = length(X.x1)/2)
Xpop <- cbind(X.x1, X.x2) %>% as.data.frame()
n.iters <- vector(length = n.sim)
all.betas <- all.naive.betas <-  matrix(nr = n.sim, nc = 12)
all.thetas <- matrix(nr = n.sim, nc = 2)
all.causal.effects <- matrix(nr = n.sim, nc = 8)
all.H001 <- all.H002 <- all.H012 <- all.H101 <- all.H102 <- all.H112 <- matrix(nr = n.sim, nc = length(H.times))
all.SEs <- all.CIs.L <- all.CIs.U <- all.IQRs <- matrix(nr = n.sim, nc = ncol(all.thetas) +
                                                        2*ncol(all.betas) + 6*length(H.times) + 8)
all.out <- vector(length = n.sim)

for(j in 1:n.sim)
  {
  a <- proc.time()
  Daniel::CatIndex(j)
  sim.df <- SimDataWeibFrail(n.sample = n.sample, params = params,
                           no.protected = no.protected, no.large = no.large,
                           cens.exp.rate = cens.rate, cens.admin = 100, round.times = F)
  my.data <- data.frame(T1 = sim.df$T1, T2 = sim.df$T2, delta1 = sim.df$delta1,
                      delta2 = sim.df$delta2, A = sim.df$A, X = sim.df$X)
  Xnames <- c("X.x1", "X.x2")
  res <- EMcausalSC(data = my.data, Xnames = Xnames, max.iter = max.iter, init.thetas = c(params$theta, params$theta))
  n.iters[j] <- res$iter
  all.betas[j, ] <- res$betas
  all.naive.betas[j, ] <- res$naive.betas
  all.thetas[j, ] <- res$thetas
  all.H001[j, ] <- res$H.step.funcs$step.A0T1(H.times)
  all.H002[j, ] <- res$H.step.funcs$step.A0T2(H.times)
  all.H012[j, ] <- res$H.step.funcs$step.A0T12(H.times)
  all.H101[j, ] <- res$H.step.funcs$step.A1T1(H.times)
  all.H102[j, ] <- res$H.step.funcs$step.A1T2(H.times)
  all.H112[j, ] <- res$H.step.funcs$step.A1T12(H.times)
  all.causal.effects[j, ] <- CalcRMST(rho = params$rho, tau = tau, n.gamma.vals = n.gamma.vals, n.sample.pers = n.sample.pers,
                                               population = Xpop, data = my.data, Xnames = Xnames, res = res, list.out = F)
  boot.obj <- boot(data = my.data, statistic = WrapEMCausalSCBoot, R = R,
                  parallel = "snow", ncpus = 10,
                  Xnames = Xnames, population = Xpop,
                  max.iter = max.iter,
                  tau = tau, rhos = params$rho,
                  n.sample.pers = n.sample.pers,
                  n.gamma.vals = n.gamma.vals,
                  init.thetas = res$thetas, H.times = H.times)
  all.out[j] <- sum(is.na(boot.obj$t[, 25]))
  all.SEs[j, ] <- apply(boot.obj$t, 2, sd, na.rm = T)
  all.IQRs[j, ] <- 20*apply(boot.obj$t, 2, IQR, na.rm = T)/27
  all.CIs <- apply(boot.obj$t, 2, quantile, c(0.025, 0.975), na.rm = T)
  all.CIs.L[j, ] <- all.CIs[1, ]
  all.CIs.U[j, ] <- all.CIs[2, ]
  b <- proc.time()
  print(b - a)
}

packageVersion("CausalSemiComp")
setwd("/home/dn84/CausalSemiComp/Results")
save.image(paste0("SP2000",scen,lett,"Cens25.RData"))
