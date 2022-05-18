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

lett <- "i"
scen <- 20001
set.seed(scen + match(lett,letters))

#####################################################################################
params <- GetScenarioParams(scenario.num = scen)
no.protected <- T
no.large <- T
cens.rate <- 0.007
n.sample <- 2000
n.sim <- 100
R <- 100
all.times <- seq(0, 50, 0.2)
SE.out.times <- seq(2, 50, 2) # times in which we calculate SE for different functions
n.times <- length(all.times)

all.delta2 <- all.delta1 <- all.etaA0 <- all.etaA1 <- vector(length = n.sim)
all.S2A0 <- all.S2A1 <- all.S1A0 <- all.S1A1 <- matrix(nr = n.sim, nc = n.times)
all.etasA0T2.le.t <- all.etasA1T2.le.t <- matrix(nr = n.sim, nc = n.times)
all.S1A0T1lT2.all <- all.S1A1T1lT2.all <- matrix(nr = n.sim, nc = n.times)
all.SEs <- matrix(nr = n.sim, nc = 8*length(SE.out.times) + 2)

true.effects <- CalcTrueCausalParams(n.sample = 1000000, params = params, all.times = all.times, adjusted = F,
                                     no.large = no.large, no.protected = no.protected)

for (j in 1:n.sim)
{
  Daniel::CatIndex(j)
  a <- proc.time()
  sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.protected = no.protected, cens.exp.rate = cens.rate)
  all.delta1[j] <- mean(sim.df$delta1)
  all.delta2[j] <- mean(sim.df$delta2)
  res <- CausalSC(T1 = sim.df$T1, T2 = sim.df$T2, delta1 = sim.df$delta1, delta2 = sim.df$delta2,
                  A = sim.df$A, all.times = all.times)
  my.data <- data.frame(T1 = sim.df$T1, T2 = sim.df$T2, delta1 = sim.df$delta1,
                        delta2 = sim.df$delta2, A = sim.df$A)
  boot.obj <- boot(data = my.data, statistic = WrapCausalSCboot, R = R,
                   parallel = "snow", ncpus = 10,
                   out.times = SE.out.times, all.times = all.times, Ltime = 0)
  all.etaA0[j] <- res$etaA0
  all.etaA1[j] <- res$etaA1
  all.S2A0[j, ] <- res$S2A0.all.times
  all.S2A1[j, ] <- res$S2A1.all.times
  all.S1A0[j, ] <- res$S1A0.all.times
  all.S1A1[j, ] <- res$S1A1.all.times
  all.etasA0T2.le.t[j, ] <- res$etasA0T2.le.t
  all.etasA1T2.le.t[j, ] <- res$etasA1T2.le.t
  all.S1A0T1lT2.all[j, ] <- res$S1A0T1lT2.all.times
  all.S1A1T1lT2.all[j, ] <- res$S1A1T1lT2.all.times
  all.SEs[j, ] <- apply(boot.obj$t, 2, sd)
  b <- proc.time()
  print(b - a)
}

packageVersion("CausalSemiComp")
save.image(paste0("NP2000",scen,"cens10",lett,".RData"))
