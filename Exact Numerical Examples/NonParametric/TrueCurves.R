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
F1.a0.ad <- causal.params$F1.a0.ad$F1.a0.ad
F1.a1.ad <- causal.params$F1.a1.ad$F1.a1.ad
F2.a0.ad <- causal.params$F2.a0.ad$F2.a0.ad
F2.a1.ad <- causal.params$F2.a1.ad$F2.a1.ad
F2.a0.nd <- causal.params$F2.a0.nd$F2.a0.nd
F2.a1.nd <- causal.params$F2.a1.nd$F2.a1.nd

###### Create df #########

df.scen1 <- data.frame(a = rep(rep(c("a=0", "a=1"), each = length(all.times)), 3),
                     Type = rep(c("F1 ad", "F2 ad", "F2 nd"), each = 2 * length(all.times)),
                     Time = rep(all.times, 6),
                     CDF = c(F1.a0.ad, F1.a1.ad, F2.a0.ad, F2.a1.ad, F2.a0.nd, F2.a1.nd))
df.scen1$Scen <- "Scenario~(I)"

rm(list = setdiff(ls(),list("df.scen1")))

set.seed(314)

all.times <- seq(0, 50, 0.2)
n.times <- length(all.times)
params <- GetScenarioParams(scenario.num = 20001)
no.protected <- T
no.large <- T
causal.params <- CalcTrueCausalParams(n.sample = 10000000,  params, all.times = all.times,
                                      no.protected = no.protected, no.large = no.large, adjusted = T)

#### True Values ####
F1.a0.ad <- causal.params$F1.a0.ad$F1.a0.ad
F1.a1.ad <- causal.params$F1.a1.ad$F1.a1.ad
F2.a0.ad <- causal.params$F2.a0.ad$F2.a0.ad
F2.a1.ad <- causal.params$F2.a1.ad$F2.a1.ad
F2.a0.nd <- causal.params$F2.a0.nd$F2.a0.nd
F2.a1.nd <- causal.params$F2.a1.nd$F2.a1.nd

###### Create df #########

df.scen2 <- data.frame(a = rep(rep(c("a=0", "a=1"), each = length(all.times)), 3),
                       Type = rep(c("F1 ad", "F2 ad", "F2 nd"), each = 2 * length(all.times)),
                       Time = rep(all.times, 6),
                       CDF = c(F1.a0.ad, F1.a1.ad, F2.a0.ad, F2.a1.ad, F2.a0.nd, F2.a1.nd))
df.scen2$Scen <- "Scenario~(II)"

rm(list = setdiff(ls(),list("df.scen1", "df.scen2")))





set.seed(314)

all.times <- seq(0, 50, 0.2)
n.times <- length(all.times)
params <- GetScenarioParams(scenario.num = 30001)
no.protected <- T
no.large <- T
causal.params <- CalcTrueCausalParams(n.sample = 10000000,  params, all.times = all.times,
                                      no.protected = no.protected, no.large = no.large, adjusted = T)

#### True Values ####
F1.a0.ad <- causal.params$F1.a0.ad$F1.a0.ad
F1.a1.ad <- causal.params$F1.a1.ad$F1.a1.ad
F2.a0.ad <- causal.params$F2.a0.ad$F2.a0.ad
F2.a1.ad <- causal.params$F2.a1.ad$F2.a1.ad
F2.a0.nd <- causal.params$F2.a0.nd$F2.a0.nd
F2.a1.nd <- causal.params$F2.a1.nd$F2.a1.nd

###### Create df #########

df.scen3 <- data.frame(a = rep(rep(c("a=0", "a=1"), each = length(all.times)), 3),
                       Type = rep(c("F1 ad", "F2 ad", "F2 nd"), each = 2 * length(all.times)),
                       Time = rep(all.times, 6),
                       CDF = c(F1.a0.ad, F1.a1.ad, F2.a0.ad, F2.a1.ad, F2.a0.nd, F2.a1.nd))
df.scen3$Scen <- "Scenario~(III)"

rm(list = setdiff(ls(),list("df.scen1", "df.scen2", "df.scen3")))


df.all <- rbind(df.scen1, df.scen2, df.scen3)

packageVersion("CausalSemiComp")
#setwd("/home/dn84/CausalSemiComp/Results")
# Save the results to play with them later, figures are created in TrueCurvesPlots.R
save.image(paste0("TrueCurves.RData"))
