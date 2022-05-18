### This Script create the illustrative bounds plots given in
### Figure 1 of Nevo & Gorfine (2020+)

load("~/Numerical Exact Examples/NonParametric/Bounds.RData")

df.all.diff.scen2$U.bound.adj <- NA
df.all.diff.scen2$L.bound.adj <- NA
df.all.diff.scen3$U.bound.adj <- NA
df.all.diff.scen3$L.bound.adj <- NA

df.all.diff.scen1$Scen <- "Scenario~(I)"
df.all.diff.scen2$Scen <- "Scenario~(II)"
df.all.diff.scen3$Scen <- "Scenario~(III)"

df.all.diff <- rbind(df.all.diff.scen1, df.all.diff.scen2, df.all.diff.scen3)


df.all.diff$U.bound.RS <- NA
df.all.diff$U.bound.RS[df.all.diff$Type=="F1 ad"] <- df.all.diff$U.bound[df.all.diff$Type=="F1 ad"]

df.all.diff$Type2[df.all.diff$Type=="F1 ad"] <- "ad:~T[1]"
df.all.diff$Type2[df.all.diff$Type=="F2 ad"] <- "ad:~T[2]"
df.all.diff$Type2[df.all.diff$Type=="F2 nd"] <- "nd:~T[2]"



ggplot(df.all.diff, aes(x = Time, y = Diff, ymin = L.bound, ymax = U.bound)) + theme_bw() +
  facet_grid(Scen~Type2, labeller = "label_parsed") +
  geom_line(size = 3, color = "black") + geom_ribbon(alpha = 0.3, fill = "red") +
  ylab("CDF difference") +
  ylim(c(-1 , 1)) +
  theme(title = element_text(size=  32),
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 22)) + # uncomment for adjusted bounds+
  geom_line(aes(y = L.bound.RS), size = 2, color = "darkgreen", linetype = "dashed") +
  geom_line(aes(y = U.bound.RS), size = 2, color = "darkgreen", linetype = "dashed") +
  geom_line(aes(x = Time, y = U.bound.adj), size = 1.5, color = "blue", lty = "11") +
  geom_line(aes(x = Time, y = L.bound.adj), size = 1.5, color = "blue", linetype = "11")


