####################################################################################
####################################################################################
# Figure 2 of Nevo & Gorfine (2020+)
####################################################################################rm(list = ls())

library(tidyverse)
library(scales)

setwd("~/CausalSemiCompReproduce/Numerical Exact Examples/SemiParametric/")
all.true.effects <- read.csv("TrueEffectsSemiParam")
n.effects <- nrow(all.true.effects)
 colnames(all.true.effects) <- c("rho", "theta", "prop.ad", "prop.nd.tau", "ATE.T2.ad",
                                 "ATE.T2.nd", "ATE.T1.ad", "med.ATE.T2.ad", "med.ATE.T2.nd",
                                 "med.ATE.T1.ad")

new.df <- data.frame(rho = rep(all.true.effects[, 1], 8), theta = (rep(all.true.effects[, 2], 8)))
new.df$Effect <- c(as.matrix(all.true.effects[, -c(1:2)]))
new.df$Stratum <- c(rep("ad",n.effects), rep("nd",n.effects), rep("ad",n.effects), rep("nd",n.effects),
                    rep("ad",2 * n.effects), rep("nd",n.effects), rep("ad",n.effects))
new.df$Type <- c(rep("Proportion", 2 * n.effects), rep("Mean", 3* n.effects), rep("Median", 3* n.effects))
new.df$Outcome <- c(rep("Strata", 2 * n.effects), rep(c(rep("T[2]", 2* n.effects), rep("T[1]", n.effects)), 2))
new.df$theta.fact <- factor(new.df$theta, levels = c(2/3, 1, 2), labels = c("theta==frac(2,3)",
                                                                            "theta==1",
                                                                            "theta==2"))
new.df$Combine <- str_c(as.character(new.df$Stratum), as.character(new.df$Outcome),sep = ":")
new.df$StratumProb <- str_c("pi", as.character(new.df$Stratum),sep = "[")
new.df$StratumProb <- str_c(new.df$StratumProb, "]",sep = "")
new.df %>% filter(Type!="Proportion") %>%
  ggplot(aes(x = rho, y = Effect, col = theta.fact, linetype = Type)) + theme_bw()  +
  facet_wrap(~Combine, nrow = 1, labeller = label_parsed) + geom_line(size = 2) +
  scale_colour_discrete(labels = parse_format()) +
  guides(col = guide_legend(title="Frailty Variance"),
         linetype = guide_legend(title = "Effect Type")) +
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        panel.spacing = unit(1.5, "lines")) +
  #legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(-7.5, 0), breaks = seq(0, -10, -2)) + xlab(expression(rho)) +
  ylab("Causal RMST Effect")

new.df %>% filter(Type=="Proportion") %>%
  ggplot(aes(x = rho, y = Effect, col = theta.fact)) + theme_bw()  +
  facet_wrap(~StratumProb, nrow = 1, labeller = label_parsed) + geom_line(size = 2) +
  scale_colour_discrete(labels = parse_format()) +
  guides(col=guide_legend(title="Frailty Variance")) +
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        panel.spacing = unit(1.5, "lines")) +
        #legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  #scale_y_continuous(breaks = seq(0, -15, -5)) +
    xlab(expression(rho)) +
  ylab("Strata proportion") + scale_y_continuous(limits = c(0,0.65))

new.df %>% dplyr::filter(Outcome == "Strata") %>% group_by(Stratum) %>%
  summarize(mn = min(Effect), mx = max(Effect)) %>% round(2)
######### Try in 1 Figure (didn't work so well) #######

new.df$title1 <- ifelse(new.df$Effect=="Proportion", "Strata~Proportion", "Causal~RMST~Effect")
new.df$title2 <- ifelse(new.df$Effect=="Proportion", new.df$StratumProb, new.df$Combine)

new.df %>%
  ggplot(aes(x = rho, y = Effect, col = theta.fact, linetype = Type)) + theme_bw()  +
  facet_wrap(title1~title2, nrow = 1, labeller = label_parsed) + geom_line(size = 2) +
  scale_colour_discrete(labels = parse_format()) +
  guides(col=guide_legend(title="Frailty Variance"),
         linetype = guide_legend(title = element_blank())) +
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, -15, -5)) + xlab(expression(rho)) +
  ylab("Causal RMST")
