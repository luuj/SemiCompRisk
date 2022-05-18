####################################################################################
####################################################################################
# Large sample bounds - choose scenario in line 17. Scenarios are translated
# to parameters in the function GetScenarioParams of the CausalSemiComp package
####################################################################################
rm(list = ls())

library(dplyr)
library(ggplot2)
library(CausalSemiComp)
####################################################################################
set.seed(314)
####################################################################################

all.times <- seq(0, 50, 0.2)
n.times <- length(all.times)
params <- GetScenarioParams(scenario.num = 10752)
no.protected <- T
no.large <- T
causal.params <- CalcTrueCausalParams(n.sample = 10000000,  params, all.times = all.times,
                                      no.protected = no.protected, no.large = no.large, adjusted = T)

#### True Values ####
F1.ad.real.diff <- causal.params$F1.a1.ad$F1.a1.ad - causal.params$F1.a0.ad$F1.a0.ad
F2.ad.real.diff <- causal.params$F2.a1.ad$F2.a1.ad - causal.params$F2.a0.ad$F2.a0.ad
F2.nd.real.diff <- causal.params$F2.a1.nd$F2.a1.nd - causal.params$F2.a0.nd$F2.a0.nd
#### Unadjusted Bounds ####
etaA0 <- causal.params$true.eta0$true.eta0
etaA1 <- causal.params$true.eta1$true.eta1
S2A0.all.times <- causal.params$true.S2A0$true.S2A0
S2A1.all.times <- causal.params$true.S2A1$true.S2A1
S1A0.all.times <- causal.params$true.S1A0$true.S1A0
S1A1.all.times <- causal.params$true.S1A1$true.S1A1
F2A0T1lT2.all.times <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2
F2A1T1lT2.all.times <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2
etasA0T2.le.t <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0
etasA1T2.le.t <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1
F1A0T1lT2.all.times <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2
F1A1T1lT2.all.times <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2
F2A0T1gT2 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2
F2A1T1gT2 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2

###### Calculate Bounds #########
ad.T2.L <- pmax(0, 1 - S2A1.all.times/etaA0) - F2A0T1lT2.all.times
ad.T2.U <- pmin(1, (1 - S2A1.all.times) * etasA1T2.le.t/etaA0) - F2A0T1lT2.all.times
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L <- pmax(0, 1 - S1A1.all.times/etaA0) - F1A0T1lT2.all.times
ad.T1.L.RS <- F1A1T1lT2.all.times - F1A0T1lT2.all.times
ad.T1.U <- pmin(1, F1A1T1lT2.all.times*etaA1/etaA0) - F1A0T1lT2.all.times


#### Adjusted Bounds ####
etaA0Z0 <- causal.params$true.eta0$true.eta0Z0
etaA0Z1 <- causal.params$true.eta0$true.eta0Z1
etaA1Z0 <- causal.params$true.eta1$true.eta1Z0
etaA1Z1 <- causal.params$true.eta1$true.eta1Z1
S2A0.all.timesZ0 <- causal.params$true.S2A0$true.S2A0Z0
S2A1.all.timesZ0 <- causal.params$true.S2A1$true.S2A1Z0
S2A0.all.timesZ1 <- causal.params$true.S2A0$true.S2A0Z1
S2A1.all.timesZ1 <- causal.params$true.S2A1$true.S2A1Z1
S1A0.all.timesZ0 <- causal.params$true.S1A0$true.S1A0Z0
S1A1.all.timesZ0 <- causal.params$true.S1A1$true.S1A1Z0
S1A0.all.timesZ1 <- causal.params$true.S1A0$true.S1A0Z1
S1A1.all.timesZ1 <- causal.params$true.S1A1$true.S1A1Z1
F2A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z0
F2A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z0
F2A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z1
F2A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z1
etasA0T2.le.tZ0 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z0
etasA1T2.le.tZ0 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z0
etasA0T2.le.tZ1 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z1
etasA1T2.le.tZ1 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z1
F1A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z0
F1A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z0
F1A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z1
F1A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z1
F2A0T1gT2Z0 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z0
F2A1T1gT2Z0 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z0
F2A0T1gT2Z1 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z1
F2A1T1gT2Z1 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z1
Zprob.ad <- causal.params$Zprob.ad # This is the probability Pr(Z = 1| T_1(0)\le T_2(0))
Zprob.nd <- causal.params$Zprob.nd # This is the probability Pr(Z = 1| T_1(1) > T_2(1))
F2A0T1lT2.all.timesZ0 <- (1 - S2A0.all.timesZ0)*etasA0T2.le.tZ0/etaA0Z0
F2A0T1lT2.all.timesZ1 <- (1 - S2A0.all.timesZ1)*etasA0T2.le.tZ1/etaA0Z1
ad.T2.L.Z0 <- pmax(0, 1 - S2A1.all.timesZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.U.Z0 <- pmin(1, (1 - S2A1.all.timesZ0) * etasA1T2.le.tZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.L.Z1 <- pmax(0, 1 - S2A1.all.timesZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
ad.T2.U.Z1 <- pmin(1, (1 - S2A1.all.timesZ1) * etasA1T2.le.tZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
nd.T2.L.Z0 <- F2A1T1gT2Z0 - pmin(1, (1 - S2A0.all.timesZ0) * (1 - etasA0T2.le.tZ0) / (1 - etaA1Z0))
nd.T2.U.Z0 <- F2A1T1gT2Z0 - pmax(0, 1  - S2A0.all.timesZ0 / (1 - etaA1Z0))
nd.T2.L.Z1 <- F2A1T1gT2Z1 - pmin(1, (1 - S2A0.all.timesZ1) * (1 - etasA0T2.le.tZ1) / (1 - etaA1Z1))
nd.T2.U.Z1 <- F2A1T1gT2Z1 - pmax(0, 1  - S2A0.all.timesZ1 / (1 - etaA1Z1))
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L.Z0 <- pmax(0, 1 - S1A1.all.timesZ0/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.L.Z1 <- pmax(0, 1 - S1A1.all.timesZ1/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.U.Z0 <- pmin(1, (1 - S1A1.all.timesZ0)/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.U.Z1 <- pmin(1, (1 - S1A1.all.timesZ1)/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.L.adj <- (1 - Zprob.ad) * ad.T1.L.Z0 + Zprob.ad * ad.T1.L.Z1
ad.T1.U.adj <- (1 - Zprob.ad) * ad.T1.U.Z0 + Zprob.ad * ad.T1.U.Z1
ad.T2.L.adj <- (1 - Zprob.ad) * ad.T2.L.Z0 + Zprob.ad * ad.T2.L.Z1
ad.T2.U.adj <- (1 - Zprob.ad) * ad.T2.U.Z0 + Zprob.ad * ad.T2.U.Z1
nd.T2.L.adj <- (1 - Zprob.nd) * nd.T2.L.Z0 + Zprob.nd * nd.T2.L.Z1
nd.T2.U.adj <- (1 - Zprob.nd) * nd.T2.U.Z0 + Zprob.nd * nd.T2.U.Z1

###### Create df #########

df.F1.diff <- data.frame(Type = "F1 ad",Time = all.times, Diff = F1.ad.real.diff,
                         U.bound = ad.T1.U,  L.bound = ad.T1.L, L.bound.RS = ad.T1.L.RS,
                         U.bound.adj = ad.T1.U.adj, L.bound.adj = ad.T1.L.adj)

df.F2.ad.diff <- data.frame(Type = "F2 ad", Time = all.times, Diff = F2.ad.real.diff,
                         U.bound = ad.T2.U,  L.bound = ad.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                         U.bound.adj = ad.T2.U.adj, L.bound.adj = ad.T2.L.adj)

df.F2.nd.diff <- data.frame(Type = "F2 nd", Time = all.times, Diff = F2.nd.real.diff,
                            U.bound = nd.T2.U,  L.bound = nd.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                            U.bound.adj = nd.T2.L.adj, L.bound.adj = nd.T2.U.adj)

df.all.diff.scen1 <- rbind(df.F1.diff, df.F2.ad.diff, df.F2.nd.diff)
df.all.diff.scen1$Scen <- 1

rm(list = setdiff(ls(),"df.all.diff.scen1"))

set.seed(314)

all.times <- seq(0, 50, 0.2)
n.times <- length(all.times)
params <- GetScenarioParams(scenario.num = 20001)
no.protected <- T
no.large <- T
causal.params <- CalcTrueCausalParams(n.sample = 10000000,  params, all.times = all.times,
                                      no.protected = no.protected, no.large = no.large, adjusted = T)

#### True Values ####
F1.ad.real.diff <- causal.params$F1.a1.ad$F1.a1.ad - causal.params$F1.a0.ad$F1.a0.ad
F2.ad.real.diff <- causal.params$F2.a1.ad$F2.a1.ad - causal.params$F2.a0.ad$F2.a0.ad
F2.nd.real.diff <- causal.params$F2.a1.nd$F2.a1.nd - causal.params$F2.a0.nd$F2.a0.nd
#### Unadjusted Bounds ####
etaA0 <- causal.params$true.eta0$true.eta0
etaA1 <- causal.params$true.eta1$true.eta1
S2A0.all.times <- causal.params$true.S2A0$true.S2A0
S2A1.all.times <- causal.params$true.S2A1$true.S2A1
S1A0.all.times <- causal.params$true.S1A0$true.S1A0
S1A1.all.times <- causal.params$true.S1A1$true.S1A1
F2A0T1lT2.all.times <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2
F2A1T1lT2.all.times <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2
etasA0T2.le.t <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0
etasA1T2.le.t <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1
F1A0T1lT2.all.times <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2
F1A1T1lT2.all.times <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2
F2A0T1gT2 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2
F2A1T1gT2 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2

###### Calculate Bounds #########
ad.T2.L <- pmax(0, 1 - S2A1.all.times/etaA0) - F2A0T1lT2.all.times
ad.T2.U <- pmin(1, (1 - S2A1.all.times) * etasA1T2.le.t/etaA0) - F2A0T1lT2.all.times
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L <- pmax(0, 1 - S1A1.all.times/etaA0) - F1A0T1lT2.all.times
ad.T1.L.RS <- F1A1T1lT2.all.times - F1A0T1lT2.all.times
ad.T1.U <- pmin(1, F1A1T1lT2.all.times*etaA1/etaA0) - F1A0T1lT2.all.times


#### Adjusted Bounds ####
etaA0Z0 <- causal.params$true.eta0$true.eta0Z0
etaA0Z1 <- causal.params$true.eta0$true.eta0Z1
etaA1Z0 <- causal.params$true.eta1$true.eta1Z0
etaA1Z1 <- causal.params$true.eta1$true.eta1Z1
S2A0.all.timesZ0 <- causal.params$true.S2A0$true.S2A0Z0
S2A1.all.timesZ0 <- causal.params$true.S2A1$true.S2A1Z0
S2A0.all.timesZ1 <- causal.params$true.S2A0$true.S2A0Z1
S2A1.all.timesZ1 <- causal.params$true.S2A1$true.S2A1Z1
S1A0.all.timesZ0 <- causal.params$true.S1A0$true.S1A0Z0
S1A1.all.timesZ0 <- causal.params$true.S1A1$true.S1A1Z0
S1A0.all.timesZ1 <- causal.params$true.S1A0$true.S1A0Z1
S1A1.all.timesZ1 <- causal.params$true.S1A1$true.S1A1Z1
F2A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z0
F2A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z0
F2A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z1
F2A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z1
etasA0T2.le.tZ0 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z0
etasA1T2.le.tZ0 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z0
etasA0T2.le.tZ1 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z1
etasA1T2.le.tZ1 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z1
F1A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z0
F1A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z0
F1A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z1
F1A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z1
F2A0T1gT2Z0 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z0
F2A1T1gT2Z0 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z0
F2A0T1gT2Z1 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z1
F2A1T1gT2Z1 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z1
Zprob.ad <- causal.params$Zprob.ad # This is the probability Pr(Z = 1| T_1(0)\le T_2(0))
Zprob.nd <- causal.params$Zprob.nd # This is the probability Pr(Z = 1| T_1(1) > T_2(1))
F2A0T1lT2.all.timesZ0 <- (1 - S2A0.all.timesZ0)*etasA0T2.le.tZ0/etaA0Z0
F2A0T1lT2.all.timesZ1 <- (1 - S2A0.all.timesZ1)*etasA0T2.le.tZ1/etaA0Z1
ad.T2.L.Z0 <- pmax(0, 1 - S2A1.all.timesZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.U.Z0 <- pmin(1, (1 - S2A1.all.timesZ0) * etasA1T2.le.tZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.L.Z1 <- pmax(0, 1 - S2A1.all.timesZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
ad.T2.U.Z1 <- pmin(1, (1 - S2A1.all.timesZ1) * etasA1T2.le.tZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
nd.T2.L.Z0 <- F2A1T1gT2Z0 - pmin(1, (1 - S2A0.all.timesZ0) * (1 - etasA0T2.le.tZ0) / (1 - etaA1Z0))
nd.T2.U.Z0 <- F2A1T1gT2Z0 - pmax(0, 1  - S2A0.all.timesZ0 / (1 - etaA1Z0))
nd.T2.L.Z1 <- F2A1T1gT2Z1 - pmin(1, (1 - S2A0.all.timesZ1) * (1 - etasA0T2.le.tZ1) / (1 - etaA1Z1))
nd.T2.U.Z1 <- F2A1T1gT2Z1 - pmax(0, 1  - S2A0.all.timesZ1 / (1 - etaA1Z1))
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L.Z0 <- pmax(0, 1 - S1A1.all.timesZ0/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.L.Z1 <- pmax(0, 1 - S1A1.all.timesZ1/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.U.Z0 <- pmin(1, (1 - S1A1.all.timesZ0)/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.U.Z1 <- pmin(1, (1 - S1A1.all.timesZ1)/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.L.adj <- (1 - Zprob.ad) * ad.T1.L.Z0 + Zprob.ad * ad.T1.L.Z1
ad.T1.U.adj <- (1 - Zprob.ad) * ad.T1.U.Z0 + Zprob.ad * ad.T1.U.Z1
ad.T2.L.adj <- (1 - Zprob.ad) * ad.T2.L.Z0 + Zprob.ad * ad.T2.L.Z1
ad.T2.U.adj <- (1 - Zprob.ad) * ad.T2.U.Z0 + Zprob.ad * ad.T2.U.Z1
nd.T2.L.adj <- (1 - Zprob.nd) * nd.T2.L.Z0 + Zprob.nd * nd.T2.L.Z1
nd.T2.U.adj <- (1 - Zprob.nd) * nd.T2.U.Z0 + Zprob.nd * nd.T2.U.Z1

###### Create df #########

df.F1.diff <- data.frame(Type = "F1 ad",Time = all.times, Diff = F1.ad.real.diff,
                         U.bound = ad.T1.U,  L.bound = ad.T1.L, L.bound.RS = ad.T1.L.RS,
                         U.bound.adj = ad.T1.U.adj, L.bound.adj = ad.T1.L.adj)

df.F2.ad.diff <- data.frame(Type = "F2 ad", Time = all.times, Diff = F2.ad.real.diff,
                            U.bound = ad.T2.U,  L.bound = ad.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                            U.bound.adj = ad.T2.U.adj, L.bound.adj = ad.T2.L.adj)

df.F2.nd.diff <- data.frame(Type = "F2 nd", Time = all.times, Diff = F2.nd.real.diff,
                            U.bound = nd.T2.U,  L.bound = nd.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                            U.bound.adj = nd.T2.L.adj, L.bound.adj = nd.T2.U.adj)

df.all.diff.scen2 <- rbind(df.F1.diff, df.F2.ad.diff, df.F2.nd.diff)
df.all.diff.scen2$Scen <- 2

rm(list = setdiff(ls(),list("df.all.diff.scen1", "df.all.diff.scen2")))


set.seed(314)

all.times <- seq(0, 50, 0.2)
n.times <- length(all.times)
params <- GetScenarioParams(scenario.num = 30001)
no.protected <- T
no.large <- T
causal.params <- CalcTrueCausalParams(n.sample = 10000000,  params, all.times = all.times,
                                      no.protected = no.protected, no.large = no.large, adjusted = T)

#### True Values ####
F1.ad.real.diff <- causal.params$F1.a1.ad$F1.a1.ad - causal.params$F1.a0.ad$F1.a0.ad
F2.ad.real.diff <- causal.params$F2.a1.ad$F2.a1.ad - causal.params$F2.a0.ad$F2.a0.ad
F2.nd.real.diff <- causal.params$F2.a1.nd$F2.a1.nd - causal.params$F2.a0.nd$F2.a0.nd
#### Unadjusted Bounds ####
etaA0 <- causal.params$true.eta0$true.eta0
etaA1 <- causal.params$true.eta1$true.eta1
S2A0.all.times <- causal.params$true.S2A0$true.S2A0
S2A1.all.times <- causal.params$true.S2A1$true.S2A1
S1A0.all.times <- causal.params$true.S1A0$true.S1A0
S1A1.all.times <- causal.params$true.S1A1$true.S1A1
F2A0T1lT2.all.times <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2
F2A1T1lT2.all.times <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2
etasA0T2.le.t <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0
etasA1T2.le.t <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1
F1A0T1lT2.all.times <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2
F1A1T1lT2.all.times <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2
F2A0T1gT2 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2
F2A1T1gT2 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2

###### Calculate Bounds #########
ad.T2.L <- pmax(0, 1 - S2A1.all.times/etaA0) - F2A0T1lT2.all.times
ad.T2.U <- pmin(1, (1 - S2A1.all.times) * etasA1T2.le.t/etaA0) - F2A0T1lT2.all.times
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L <- pmax(0, 1 - S1A1.all.times/etaA0) - F1A0T1lT2.all.times
ad.T1.L.RS <- F1A1T1lT2.all.times - F1A0T1lT2.all.times
ad.T1.U <- pmin(1, F1A1T1lT2.all.times*etaA1/etaA0) - F1A0T1lT2.all.times


#### Adjusted Bounds ####
etaA0Z0 <- causal.params$true.eta0$true.eta0Z0
etaA0Z1 <- causal.params$true.eta0$true.eta0Z1
etaA1Z0 <- causal.params$true.eta1$true.eta1Z0
etaA1Z1 <- causal.params$true.eta1$true.eta1Z1
S2A0.all.timesZ0 <- causal.params$true.S2A0$true.S2A0Z0
S2A1.all.timesZ0 <- causal.params$true.S2A1$true.S2A1Z0
S2A0.all.timesZ1 <- causal.params$true.S2A0$true.S2A0Z1
S2A1.all.timesZ1 <- causal.params$true.S2A1$true.S2A1Z1
S1A0.all.timesZ0 <- causal.params$true.S1A0$true.S1A0Z0
S1A1.all.timesZ0 <- causal.params$true.S1A1$true.S1A1Z0
S1A0.all.timesZ1 <- causal.params$true.S1A0$true.S1A0Z1
S1A1.all.timesZ1 <- causal.params$true.S1A1$true.S1A1Z1
F2A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z0
F2A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z0
F2A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A0T1leT2$true.S2A0T1leT2Z1
F2A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S2A1T1leT2$true.S2A1T1leT2Z1
etasA0T2.le.tZ0 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z0
etasA1T2.le.tZ0 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z0
etasA0T2.le.tZ1 <- causal.params$true.etas.T2leT.0$true.etas.T2leT.0.Z1
etasA1T2.le.tZ1 <- causal.params$true.etas.T2leT.1$true.etas.T2leT.1.Z1
F1A0T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z0
F1A1T1lT2.all.timesZ0 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z0
F1A0T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A0T1leT2$true.S1A0T1leT2Z1
F1A1T1lT2.all.timesZ1 <- 1 - causal.params$true.S1A1T1leT2$true.S1A1T1leT2Z1
F2A0T1gT2Z0 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z0
F2A1T1gT2Z0 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z0
F2A0T1gT2Z1 <- 1 - causal.params$true.S2A0T1gT2$true.S2A0T1gT2Z1
F2A1T1gT2Z1 <- 1 - causal.params$true.S2A1T1gT2$true.S2A1T1gT2Z1
Zprob.ad <- causal.params$Zprob.ad # This is the probability Pr(Z = 1| T_1(0)\le T_2(0))
Zprob.nd <- causal.params$Zprob.nd # This is the probability Pr(Z = 1| T_1(1) > T_2(1))
F2A0T1lT2.all.timesZ0 <- (1 - S2A0.all.timesZ0)*etasA0T2.le.tZ0/etaA0Z0
F2A0T1lT2.all.timesZ1 <- (1 - S2A0.all.timesZ1)*etasA0T2.le.tZ1/etaA0Z1
ad.T2.L.Z0 <- pmax(0, 1 - S2A1.all.timesZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.U.Z0 <- pmin(1, (1 - S2A1.all.timesZ0) * etasA1T2.le.tZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
ad.T2.L.Z1 <- pmax(0, 1 - S2A1.all.timesZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
ad.T2.U.Z1 <- pmin(1, (1 - S2A1.all.timesZ1) * etasA1T2.le.tZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
nd.T2.L.Z0 <- F2A1T1gT2Z0 - pmin(1, (1 - S2A0.all.timesZ0) * (1 - etasA0T2.le.tZ0) / (1 - etaA1Z0))
nd.T2.U.Z0 <- F2A1T1gT2Z0 - pmax(0, 1  - S2A0.all.timesZ0 / (1 - etaA1Z0))
nd.T2.L.Z1 <- F2A1T1gT2Z1 - pmin(1, (1 - S2A0.all.timesZ1) * (1 - etasA0T2.le.tZ1) / (1 - etaA1Z1))
nd.T2.U.Z1 <- F2A1T1gT2Z1 - pmax(0, 1  - S2A0.all.timesZ1 / (1 - etaA1Z1))
nd.T2.L <- F2A1T1gT2  - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
nd.T2.U <- F2A1T1gT2 -  pmax(0, 1 - S2A0.all.times / (1 - etaA1))
ad.T1.L.Z0 <- pmax(0, 1 - S1A1.all.timesZ0/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.L.Z1 <- pmax(0, 1 - S1A1.all.timesZ1/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.U.Z0 <- pmin(1, (1 - S1A1.all.timesZ0)/etaA0Z0) - F1A0T1lT2.all.timesZ0
ad.T1.U.Z1 <- pmin(1, (1 - S1A1.all.timesZ1)/etaA0Z1) - F1A0T1lT2.all.timesZ1
ad.T1.L.adj <- (1 - Zprob.ad) * ad.T1.L.Z0 + Zprob.ad * ad.T1.L.Z1
ad.T1.U.adj <- (1 - Zprob.ad) * ad.T1.U.Z0 + Zprob.ad * ad.T1.U.Z1
ad.T2.L.adj <- (1 - Zprob.ad) * ad.T2.L.Z0 + Zprob.ad * ad.T2.L.Z1
ad.T2.U.adj <- (1 - Zprob.ad) * ad.T2.U.Z0 + Zprob.ad * ad.T2.U.Z1
nd.T2.L.adj <- (1 - Zprob.nd) * nd.T2.L.Z0 + Zprob.nd * nd.T2.L.Z1
nd.T2.U.adj <- (1 - Zprob.nd) * nd.T2.U.Z0 + Zprob.nd * nd.T2.U.Z1

###### Create df #########

df.F1.diff <- data.frame(Type = "F1 ad",Time = all.times, Diff = F1.ad.real.diff,
                         U.bound = ad.T1.U,  L.bound = ad.T1.L, L.bound.RS = ad.T1.L.RS,
                         U.bound.adj = ad.T1.U.adj, L.bound.adj = ad.T1.L.adj)

df.F2.ad.diff <- data.frame(Type = "F2 ad", Time = all.times, Diff = F2.ad.real.diff,
                            U.bound = ad.T2.U,  L.bound = ad.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                            U.bound.adj = ad.T2.U.adj, L.bound.adj = ad.T2.L.adj)

df.F2.nd.diff <- data.frame(Type = "F2 nd", Time = all.times, Diff = F2.nd.real.diff,
                            U.bound = nd.T2.U,  L.bound = nd.T2.L, L.bound.RS = rep(NA, length(F2.ad.real.diff)),
                            U.bound.adj = nd.T2.L.adj, L.bound.adj = nd.T2.U.adj)

df.all.diff.scen3 <- rbind(df.F1.diff, df.F2.ad.diff, df.F2.nd.diff)
df.all.diff.scen3$Scen <- 3

rm(list = setdiff(ls(),list("df.all.diff.scen1", "df.all.diff.scen2", "df.all.diff.scen3")))


df.all.diff <- rbind(df.all.diff.scen1, df.all.diff.scen2, df.all.diff.scen3)

packageVersion("CausalSemiComp")
# Save the results to play with them later, figures are created in BoundsPlots.R
save.image(paste0("Bounds.RData"))
