####################################################################################
####################################################################################
# Calculate True Effects for different rho and theta values
####################################################################################
rm(list = ls())

library(dplyr)
library(ggplot2)
library(CausalSemiComp)
####################################################################################
set.seed(314)

# Define the target population
X.x1 <- rep(seq(-2, 2, 0.5), times = 2)
X.x2 <- rep(0:1, each = length(X.x1)/2)
Xpop <- cbind(X.x1, X.x2) %>% as.data.frame()

# Define tau
tau <- 15

# All scenarios considered
all.scens.5xxx1 <- seq(50001, 51001, 100)
all.scens.5xxx2 <- seq(50002, 51002, 100)
all.scens.5xxx3 <- seq(50003, 51003, 100)
all.scens <- c(all.scens.5xxx1, all.scens.5xxx2, all.scens.5xxx3)
n.scens <- length(all.scens)
# Make space
all.true.effects <- matrix(nr = length(all.scens), nc = 1 + 1 + 2 + 6)


for (i in 1:length(all.scens))
  {
  a <- proc.time()
  Daniel::CatIndex(i)
  scen <- all.scens[i]
  params <- GetScenarioParams(scenario.num = scen)
  all.true.effects[i, 1] <- params$rho
  all.true.effects[i, 2] <- params$theta
  scen.true.effects <- CalcTrueCausalParams(n.sample = 10^8, params = params,
                                     no.large = F, no.protected = F, round.times = F,
                                     X = Xpop, tau = tau, RMST.only = T)
  all.true.effects[i, 3] <- scen.true.effects$prop.ad.tau
  all.true.effects[i, 4] <- scen.true.effects$prop.nd.tau
  all.true.effects[i, 5] <- scen.true.effects$ATE.T2.ad
  all.true.effects[i, 6] <- scen.true.effects$ATE.T2.nd
  all.true.effects[i, 7] <- scen.true.effects$ATE.T1.ad
  all.true.effects[i, 8] <- scen.true.effects$med.ATE.T2.ad
  all.true.effects[i, 9] <- scen.true.effects$med.ATE.T2.nd
  all.true.effects[i, 10] <- scen.true.effects$med.ATE.T1.ad
  b <- proc.time()
  gc()
  print(b - a)
}

## Save true effects in a csv table file
df.all.true.effects <- as.data.frame(all.true.effects)
colnames(df.all.true.effects) <- c("rho","theta","prop.ad","prop.nd.tau",
                      "ATE.T2.ad","ATE.T2.nd","ATE.T1.ad",
                      "med.ATE.T2.ad","med.ATE.T2.nd","med.ATE.T1.ad")

write.csv(df.all.true.effects, "TrueEffectsLargeNnoRound", row.names = F)

packageVersion("CausalSemiComp")
#setwd("/home/dn84/CausalSemiComp/Results")
setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Results/")
save.image("CausalEffects50000LargerNnoRound.RData")
