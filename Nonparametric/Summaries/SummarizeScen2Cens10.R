####################################################################################
####################################################################################
# CausalSemiComp
# Summarize nonparametric estimation
####################################################################################
rm(list = ls())

library(dplyr)
library(xtable)
library(CausalSemiComp)

setwd("CausalSemiCompReproduce/Nonparametric/Results/Scenario 2/cens10/")
my.files <- list.files()
my.files <- my.files[substr(my.files, 1,2)=="NP"]
report.times <- c(10,  20, 30)
n.report.times <- length(report.times)
## 4: Scnario, cesnoring, delta1, delta2
# n.report.times - number of times for which we report time-varying estimands
# 2 + ... - the 2 is etaA0 etaA1
# 2* (...) - the 2 is one for SE, one for est.
all.results <- matrix(nr = 100000, nc = 4 + 2*(2 + 8*n.report.times))
free <- 1
# see SE.out.times
for(ii in 1:length(my.files))
{
  load(paste0(my.files[ii]))
  scenario <- as.numeric(substr(my.files[ii], 7, 7))
  cens <- as.numeric(substr(my.files[ii], 16, 17))
  range.sims <- free:(free + n.sim - 1)
  ind.times.all <- all.times %in% report.times
  ind.times.se <- SE.out.times %in% report.times
  all.results[range.sims, 1] <- scenario
  all.results[range.sims, 2] <- cens
  all.results[range.sims, 3] <- all.delta1
  all.results[range.sims, 4] <- all.delta2
  all.results[range.sims, 5] <- all.etaA0
  all.results[range.sims, 6] <- all.etaA1
  all.results[range.sims, (6 + 1):(6 + n.report.times)] <- all.S2A0[, ind.times.all]
  all.results[range.sims, (6 + n.report.times + 1):(6 + 2 *n.report.times)] <- all.S2A1[, ind.times.all]
  all.results[range.sims, (6 + 2 * n.report.times + 1):(6 + 3 * n.report.times)] <- all.S1A0[, ind.times.all]
  all.results[range.sims, (6 + 3 * n.report.times + 1):(6 + 4 * n.report.times)] <- all.S1A1[, ind.times.all]
  all.results[range.sims, (6 + 4 * n.report.times + 1):(6 + 5 * n.report.times)] <- all.etasA0T2.le.t[, ind.times.all]
  all.results[range.sims, (6 + 5 * n.report.times + 1):(6 + 6 * n.report.times)] <- all.etasA1T2.le.t[, ind.times.all]
  all.results[range.sims, (6 + 6 * n.report.times + 1):(6 + 7 * n.report.times)] <- all.S1A0T1lT2.all[, ind.times.all]
  all.results[range.sims, (6 + 7 * n.report.times + 1):(6 + 8 * n.report.times)] <- all.S1A1T1lT2.all[, ind.times.all]
  all.results[range.sims, (6 + 8 * n.report.times + 1):ncol(all.results)] <-
    all.SEs[, c(T, T, rep(ind.times.se, 8))]
  free <- free + n.sim
}
all.results <- all.results[1:(free -1),]
nrow(all.results)

colnames(all.results) <- c("Scenario", "Censoring", "delta1", "delta2",
                           "etaA0", "etaA1",
                           paste0("S2A0", report.times), paste0("S2A1", report.times),
                           paste0("S1A0", report.times), paste0("S1A1", report.times),
                           paste0("etaA0T2.le.t", report.times), paste0("etaA1T2.le.t", report.times),
                           paste0("S1A0T1lT2", report.times), paste0("S1A1T1lT2", report.times),
                           "SE.etaA0", "SE.etaA1",
                           paste0("SE.S2A0", report.times), paste0("SE.S2A1", report.times),
                           paste0("SE.S1A0", report.times), paste0("SE.S1A1", report.times),
                           paste0("SE.etaA0T2.le.t", report.times), paste0("SE.etaA1T2.le.t", report.times),
                           paste0("SE.S1A0T1lT2", report.times), paste0("SE.S1A1T1lT2", report.times))
all.results <- as.data.frame(all.results)

causal.params <- CalcTrueCausalParams(n.sample = 10000000, all.times = all.times, no.protected = T, no.large = T,
                                      params = params) # This could be very long - decrease n.sample if it is not needed for reporting
######## Table 1 #########
all.times[ind.times.all]
true.vals1 <- c(causal.params$true.eta0, causal.params$true.eta1,
                causal.params$true.S1A0[ind.times.all][c(1,3)],
                causal.params$true.S1A1[ind.times.all][c(1,3)])
true.vals1 %>% round(3)

all.results.tab1 <- all.results %>% select(etaA0, etaA1,
                                           S1A010, S1A030, S1A110, S1A130,
                                           SE.etaA0, SE.etaA1, SE.S1A010, SE.S1A030,
                                           SE.S1A110, SE.S1A130)

all.results.tab1 <- all.results.tab1 %>%
  mutate(CI.L.etaA0 = etaA0 - 1.96*SE.etaA0, CI.H.etaA0 = etaA0 + 1.96*SE.etaA0,
  CI.L.etaA1 = etaA1 - 1.96*SE.etaA1, CI.H.etaA1 = etaA1 + 1.96*SE.etaA1,
  CI.L.S1A010 = S1A010 - 1.96*SE.S1A010, CI.H.S1A010 = S1A010 + 1.96*SE.S1A010,
  CI.L.S1A030 = S1A030 - 1.96*SE.S1A030, CI.H.S1A030 = S1A030 + 1.96*SE.S1A030,
  CI.L.S1A110 = S1A110 - 1.96*SE.S1A110, CI.H.S1A110 = S1A110 + 1.96*SE.S1A110,
  CI.L.S1A130 = S1A130 - 1.96*SE.S1A130, CI.H.S1A130 = S1A130 + 1.96*SE.S1A130)

res.sum1 <- all.results.tab1  %>% summarise(
  m.etaA0 = mean(etaA0),
  m.etaA1 = mean(etaA1),
  m.S1A010 = mean(S1A010),
  m.S1A030 = mean(S1A030),
  m.S2A110 = mean(S1A110),
  m.S2A130 = mean(S1A130),
  emp.sd.etaA0 = sd(etaA0),
  emp.sd.etaA1 = sd(etaA1),
  emp.sd.S1A010 = sd(S1A010),
  emp.sd.S1A030 = sd(S1A030),
  emp.sd.S1A110 = sd(S1A110),
  emp.sd.S1A130 = sd(S1A130),
  est.se.etaA0 = mean(SE.etaA0),
  est.se.etaA1 = mean(SE.etaA1),
  est.se.S1A010 = mean(SE.S1A010),
  est.se.S1A030 = mean(SE.S1A030),
  est.se.S1A110 = mean(SE.S1A110),
  est.se.S1A130 = mean(SE.S1A130),
  cover.etaA0 = mean(CI.L.etaA0 < true.vals1[1] & CI.H.etaA0 > true.vals1[1]),
  cover.etaA1 = mean(CI.L.etaA1 < true.vals1[2] & CI.H.etaA1 > true.vals1[2]),
  cover.S1A010 = mean(CI.L.S1A010 < true.vals1[3] & CI.H.S1A010 > true.vals1[3]),
  cover.S1A030 = mean(CI.L.S1A030 < true.vals1[4] & CI.H.S1A030 > true.vals1[4]),
  cover.S1A110 = mean(CI.L.S1A110 < true.vals1[5] & CI.H.S1A110 > true.vals1[5]),
  cover.S1A130 = mean(CI.L.S1A130 < true.vals1[6] & CI.H.S1A130 > true.vals1[6]),
  CI.length.etaA0 = mean( CI.H.etaA0 - CI.L.etaA0),
  CI.length.etaA1 =  mean(CI.H.etaA1 - CI.L.etaA1),
  CI.length.S1A010 =  mean(CI.H.S1A010 - CI.L.S1A010),
  CI.length.S1A030 =  mean(CI.H.S1A030 - CI.L.S1A030),
 CI.length.S1A130 =  mean(CI.H.S1A130 - CI.L.S1A130),
 CI.length.S1A130 =  mean(CI.H.S1A130 - CI.L.S1A130),
 )


Tab1 <- rbind(true.vals1, as.numeric(select(res.sum1, starts_with("m"))),
      as.numeric(select(res.sum1, starts_with("emp"))),
      as.numeric(select(res.sum1, starts_with("est"))),
      as.numeric(select(res.sum1, starts_with("cover")))) %>% round(3)
Tab1
rownames(Tab1) <-  c("True", "Mean.EST", "EMP.SD", "EST.SE", "CI95")
Tab1 %>% xtable(digits = 3)

######## Table 2 ##########
true.vals2 <- c(causal.params$true.etas.T2leT.0[ind.times.all][c(1,2)],
                causal.params$true.etas.T2leT.1[ind.times.all][c(1,2)],
                causal.params$true.S1A0T1leT2[ind.times.all][c(1,2)],
                causal.params$true.S1A1T1leT2[ind.times.all][c(1,2)])
true.vals2 %>% round(3)

all.results.tab2 <- all.results %>% select(etaA0T2.le.t10, etaA0T2.le.t20,
                                           etaA1T2.le.t10, etaA1T2.le.t20,
                                           S1A0T1lT210, S1A0T1lT220,
                                           S1A1T1lT210, S1A1T1lT220,
                                           SE.etaA0T2.le.t10, SE.etaA0T2.le.t20,
                                           SE.etaA1T2.le.t10, SE.etaA1T2.le.t20,
                                           SE.S1A0T1lT210, SE.S1A0T1lT220,
                                           SE.S1A1T1lT210, SE.S1A1T1lT220)

all.results.tab2 <- all.results.tab2 %>%
  mutate(CI.L.etaA0T2.le.t10 = etaA0T2.le.t10 - 1.96*SE.etaA0T2.le.t10,
         CI.H.etaA0T2.le.t10 = etaA0T2.le.t10 + 1.96*SE.etaA0T2.le.t10,
         CI.L.etaA0T2.le.t20 = etaA0T2.le.t20 - 1.96*SE.etaA0T2.le.t20,
         CI.H.etaA0T2.le.t20 = etaA0T2.le.t20 + 1.96*SE.etaA0T2.le.t20,
         CI.L.etaA1T2.le.t10 = etaA1T2.le.t10 - 1.96*SE.etaA1T2.le.t10,
         CI.H.etaA1T2.le.t10 = etaA1T2.le.t10 + 1.96*SE.etaA1T2.le.t10,
         CI.L.etaA1T2.le.t20 = etaA1T2.le.t20 - 1.96*SE.etaA1T2.le.t20,
         CI.H.etaA1T2.le.t20 = etaA1T2.le.t20 + 1.96*SE.etaA1T2.le.t20,
         CI.L.S1A0T1lT210 = S1A0T1lT210 - 1.96*SE.S1A0T1lT210,
         CI.H.S1A0T1lT210 = S1A0T1lT210 + 1.96*SE.S1A0T1lT210,
         CI.L.S1A0T1lT220 = S1A0T1lT220 - 1.96*SE.S1A0T1lT220,
         CI.H.S1A0T1lT220 = S1A0T1lT220 + 1.96*SE.S1A0T1lT220,
         CI.L.S1A1T1lT210 = S1A1T1lT210 - 1.96*SE.S1A1T1lT210,
         CI.H.S1A1T1lT210 = S1A1T1lT210 + 1.96*SE.S1A1T1lT210,
         CI.L.S1A1T1lT220 = S1A1T1lT220 - 1.96*SE.S1A1T1lT220,
         CI.H.S1A1T1lT220 = S1A1T1lT220 + 1.96*SE.S1A1T1lT220)


res.sum2 <- all.results.tab2  %>% summarise(
  m.etaA0T2.le.t10 = mean(etaA0T2.le.t10),
  m.etaA0T2.le.t20 = mean(etaA0T2.le.t20),
  m.etaA1T2.le.t10 = mean(etaA1T2.le.t10),
  m.etaA1T2.le.t20 = mean(etaA1T2.le.t20),
  m.S1A0T1lT210 = mean(S1A0T1lT210),
  m.S1A0T1lT220 = mean(S1A0T1lT220),
  m.S1A1T1lT210 = mean(S1A1T1lT210),
  m.S1A1T1lT220 = mean(S1A1T1lT220),
  emp.sd.etaA0T2.le.t10 = sd(etaA0T2.le.t10),
  emp.sd.etaA0T2.le.t20 = sd(etaA0T2.le.t20),
  emp.sd.etaA1T2.le.t10 = sd(etaA1T2.le.t10),
  emp.sd.etaA1T2.le.t20 = sd(etaA1T2.le.t20),
  emp.sd.S1A0T1lT210 = sd(S1A0T1lT210),
  emp.sd.S1A0T1lT220 = sd(S1A0T1lT220),
  emp.sd.S1A1T1lT210 = sd(S1A1T1lT210),
  emp.sd.S1A1T1lT220 = sd(S1A1T1lT220),
  est.se.etaA0T2.le.t10 = mean(SE.etaA0T2.le.t10),
  est.se.etaA0T2.le.t20 = mean(SE.etaA0T2.le.t20),
  est.se.etaA1T2.le.t10 = mean(SE.etaA1T2.le.t10),
  est.se.etaA1T2.le.t20 = mean(SE.etaA1T2.le.t20),
  est.se.S1A0T1lT210 = mean(SE.S1A0T1lT210),
  est.se.S1A0T1lT220 = mean(SE.S1A0T1lT220),
  est.se.S1A1T1lT210 = mean(SE.S1A1T1lT210),
  est.se.S1A1T1lT220 = mean(SE.S1A1T1lT220),
  cover.etaA0T2.le.t10 = mean(CI.L.etaA0T2.le.t10 < true.vals2[1] & CI.H.etaA0T2.le.t10 > true.vals2[1]),
  cover.etaA0T2.le.t20 = mean(CI.L.etaA0T2.le.t20 < true.vals2[2] & CI.H.etaA0T2.le.t20 > true.vals2[2]),
  cover.etaA1T2.le.t10 = mean(CI.L.etaA1T2.le.t10 < true.vals2[3] & CI.H.etaA1T2.le.t10 > true.vals2[3]),
  cover.etaA1T2.le.t20 = mean(CI.L.etaA1T2.le.t20 < true.vals2[4] & CI.H.etaA1T2.le.t20 > true.vals2[4]),
  cover.S1A0T1lT210 = mean(CI.L.S1A0T1lT210 < true.vals2[5] & CI.H.S1A0T1lT210 > true.vals2[5]),
  cover.S1A0T1lT220 = mean(CI.L.S1A0T1lT220 < true.vals2[6] & CI.H.S1A0T1lT220 > true.vals2[6]),
  cover.S1A1T1lT210 = mean(CI.L.S1A1T1lT210 < true.vals2[7] & CI.H.S1A1T1lT210 > true.vals2[7]),
  cover.S1A1T1lT220 = mean(CI.L.S1A1T1lT220 < true.vals2[8] & CI.H.S1A1T1lT220 > true.vals2[8]))


Tab2 <- rbind(true.vals2, as.numeric(select(res.sum2, starts_with("m"))),
              as.numeric(select(res.sum2, starts_with("emp"))),
              as.numeric(select(res.sum2, starts_with("est"))),
              as.numeric(select(res.sum2, starts_with("cover")))) %>% round(3)
Tab2
rownames(Tab2) <-  c("True", "Mean.EST", "EMP.SD", "EST.SE", "CI95")
Tab2 %>% xtable(digits = 3)

