### This Script create the true curves plot given in
### Figure D.2 of the Supplementary Materials of Nevo & Gorfine (2020+)

load("~/Numerical Exact Examples/TrueCurves.RData")

df.all$Type2[df.all$Type=="F1 ad"] <- "ad:~T[1]"
df.all$Type2[df.all$Type=="F2 ad"] <- "ad:~T[2]"
df.all$Type2[df.all$Type=="F2 nd"] <- "nd:~T[2]"

ggplot(df.all, aes(x = Time, y = CDF, col = a)) + theme_bw() +
  facet_grid(Scen~Type2, labeller = "label_parsed") +
  geom_line(size = 1.75) +
  ylab("CDF") +
  ylim(c(0 , 1)) +
  theme(title = element_text(size=  32),
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 22),
        legend.text = element_text(size = 26))


