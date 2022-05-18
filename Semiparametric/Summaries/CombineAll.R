########################################################################################################
# Combine semi-parameteric simulation results into a single table
# For :"Causal Inference for Semi-competing risks data" by Nevo and Gorfine
########################################################################################################
library(dplyr)
library(xtable)
library(Daniel)
rm(list = ls())

n.iters.all <- 1000 * 3 * 2 * 3
thetas <- matrix(nr = n.iters.all, nc = 2)
betas <- matrix(nr = n.iters.all, nc = 12)
naive.betas <- matrix(nr = n.iters.all, nc = 12)
H001 <- H002 <- H012 <- H101 <- H102 <- H112 <- matrix(nr = n.iters.all, nc = 100)
causal.effects <- matrix(nr = n.iters.all, nc = 8)
CI.L <- CI.U <- SEs <- matrix(nr = n.iters.all, nc = 634)
cens <- scenario <- cens.check <- vector(length = n.iters.all)

free <- 1

all.libs <- c("cens5/50001/", "cens25/50001/", "cens5/50002/","cens25/50002/", "cens5/50003/","cens25/50003/",
              "cens5/50501/", "cens25/50501/", "cens5/50502/","cens25/50502/", "cens5/50503/","cens25/50503/",
              "cens5/51001/", "cens25/51001/", "cens5/51002/","cens25/51002/", "cens5/51003/","cens25/51003/")

n.libs <- length(all.libs)

for (jj in 1:n.libs)
{
setwd(paste0("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Results/50000/", all.libs[jj]))

all.files <- list.files()
all.files <- all.files[substr(all.files, 1, 2)=="SP"]
n.files <- length(all.files)
for (ii in 1:n.files)
{
  load(paste0(all.files[ii]))
  st <- free
  end <- free + nrow(all.thetas) - 1
  scenario[st:end] <- scen
  cens[st:end] <- ifelse(substr(all.libs[jj], 5 ,5) %>% as.numeric()==5, 5, 25)
  cens.check[st:end] <- cens.rate
  thetas[st:end, ] <- all.thetas
  betas[st:end, ] <- all.betas
  naive.betas[st:end, ] <- all.naive.betas
  causal.effects[st:end, ] <- all.causal.effects
  H001[st:end, ] <- all.H001
  H002[st:end, ] <- all.H002
  H012[st:end, ] <- all.H012
  H101[st:end, ] <- all.H101
  H102[st:end, ] <- all.H102
  H112[st:end, ] <- all.H112
  SEs[st:end, ] <- all.SEs
  CI.L[st:end, ] <- all.CIs.L
  CI.U[st:end, ] <- all.CIs.U
  free <- end + 1
}}




df.res <- cbind(scenario, cens, naive.betas, betas, thetas, H001, H002, H012, H101, H102, H112, causal.effects, SEs, CI.L, CI.U) %>% as.data.frame
#colnames(df.res) <-
beta.names <- c("beta.a0.01.1", "beta.a0.01.2", "beta.a0.02.1", "beta.a0.02.2", "beta.a0.12.1", "beta.a0.12.2",
                "beta.a1.01.1", "beta.a1.01.2", "beta.a1.02.1", "beta.a1.02.2", "beta.a1.12.1", "beta.a1.12.2")
H.names <- c(paste0("H001.t", H.times), paste0("H002.t", H.times), paste0("H012.t", H.times),
             paste0("H101.t", H.times), paste0("H102.t", H.times), paste0("H112.t", H.times))
all.param.names <-
  c(paste0("naive.", beta.names), beta.names, "theta0", "theta1", H.names,  "ATE.T2.ad", "ATE.T2.nd",
    "ATE.T1.ad", "med.ATE.T2.ad", "med.ATE.T2.nd", "med.ATE.T1.ad", "prop.ad.tau", "prop.nd.tau")
colnames(df.res) <- c("scenario", "cens", all.param.names, paste0("SE.", all.param.names), paste0("CI.L.perc.", all.param.names),
                      paste0("CI.U.perc.", all.param.names))


###### Add true causal effects ##########
setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Results/50000/")
all.true.effects <- read.csv("TrueEffectsLargeNnoRound")
colnames(all.true.effects) <- c("rho", "theta", "true.prop.ad.tau", "true.prop.nd.tau", "true.ATE.T2.ad",
                                "true.ATE.T2.nd", "true.ATE.T1.ad", "true.med.ATE.T2.ad", "true.med.ATE.T2.nd",
                                "true.med.ATE.T1.ad")

true.effects.50501 <-  all.true.effects %>% as.data.frame %>% filter(rho==0.5 & round(theta, 2)==0.67) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.50502 <-  all.true.effects %>% as.data.frame %>% filter(rho==0.5 & round(theta, 2)==1) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.50503 <-  all.true.effects %>% as.data.frame %>% filter(rho==0.5 & theta==2) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.50001 <-  all.true.effects %>% as.data.frame %>% filter(rho==0 & round(theta, 2)==0.67) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.50002 <-  all.true.effects %>% as.data.frame %>% filter(rho==0 & round(theta, 2)==1) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.50003 <-  all.true.effects %>% as.data.frame %>% filter(rho==0 & theta==2) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.51001 <-  all.true.effects %>% as.data.frame %>% filter(rho==1 & round(theta, 2)==0.67) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.51002 <-  all.true.effects %>% as.data.frame %>% filter(rho==1 & round(theta, 2)==1) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)
true.effects.51003 <-  all.true.effects %>% as.data.frame %>% filter(rho==1 & theta==2) %>%
  select(true.ATE.T1.ad, true.med.ATE.T1.ad, true.ATE.T2.ad, true.med.ATE.T2.ad, true.ATE.T2.nd, true.med.ATE.T2.nd,
         true.prop.ad.tau, true.prop.nd.tau)


true.effects.for.df <- rbind(bind_rows(replicate(2000, true.effects.50001, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.50002, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.50003, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.50501, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.50502, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.50503, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.51001, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.51002, simplify = FALSE)),
                             bind_rows(replicate(2000, true.effects.51003, simplify = FALSE)))

###### Add true parameters values ##########

### betas

true.betas <- as.data.frame(t(c(params$beta.a0.01, params$beta.a0.02, params$beta.a0.12,
                params$beta.a1.01, params$beta.a1.02, params$beta.a1.12)))
colnames(true.betas) <- paste0("true.", beta.names)
true.betas.for.df <- bind_rows(replicate(18000, true.betas, simplify = FALSE))

#### Baseline Hazards

true.H <- as.data.frame(t(c(-pweibull(H.times, shape = params$base.weib.shape.a0.01, scale = params$base.weib.scale.a0.01,
           lower.tail = F, log.p = T),
  -pweibull(H.times, shape = params$base.weib.shape.a0.02, scale = params$base.weib.scale.a0.02,
            lower.tail = F, log.p = T),
  -pweibull(H.times, shape = params$base.weib.shape.a0.12, scale = params$base.weib.scale.a0.12,
            lower.tail = F, log.p = T),
  -pweibull(H.times, shape = params$base.weib.shape.a1.01, scale = params$base.weib.scale.a1.01,
            lower.tail = F, log.p = T),
  -pweibull(H.times, shape = params$base.weib.shape.a1.02, scale = params$base.weib.scale.a1.02,
            lower.tail = F, log.p = T),
  -pweibull(H.times, shape = params$base.weib.shape.a1.12, scale = params$base.weib.scale.a1.12,
            lower.tail = F, log.p = T))))

colnames(true.H) <- paste0("true.", H.names)

true.H.for.df <- bind_rows(replicate(18000, true.H, simplify = FALSE))

###### Combine everything together

df <- cbind(arrange(df.res, scenario), true.effects.for.df, true.betas.for.df, true.H.for.df)
df$scenario %>% table
df$cens %>% table

setwd("/Users/danielnevo/Dropbox/CausalSemiComp/Sims/SemiParam/Summaries/AllResults/")
write.csv(df, "AllSimResults.csv", row.names = F)
