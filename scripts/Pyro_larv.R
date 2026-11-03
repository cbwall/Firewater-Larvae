## Pyro Larvae

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

# use pacman to load all the packages you are missing!
pacman::p_load('knitr', 'lme4', 'visreg', 'performance', 'see', "lmerTest", 'tidyverse', 'effects', 'plyr', 'dplyr', 'plotrix', 'car', "gridExtra", "cowplot", "tools", "mgcv", "gratia", "MASS", "stats", "tidymv", "sjstats", "coin", "emmeans", "ggplot2", "mda", "nortest", "reshape2", "gmm", "propagate","ggmap", "RgoogleMaps", "MixSIAR","GGally", "ggbiplot", "ggcorrplot", "vegan", "pairwiseAdonis", 'ggpubr')


df.PyLv<-read.csv("data/Pyro_larv.csv")
head(df.PyLv)

df.PyLv$level<-factor(df.PyLv$level, levels=c("control", "low", "medium", "high"))
df.PyLv$treatment<-factor(df.PyLv$treatment, levels=c("control","inorganic", "burned"))
df.PyLv$days<-factor(df.PyLv$days, levels=c("2", "5", "7"))

str(df.PyLv)

larv.plot<-ggplot(df.PyLv, aes(x=level, y=larvae))+
  geom_boxplot(fill="coral", alpha=0.9)+
  geom_jitter(color="brown", size=2, alpha=0.5) +
  facet_wrap(~days + treatment)+
  ylab("Surviving Larvae")+
  theme_classic()

sum<-aggregate(larvae~days*treatment*level, data=df.PyLv, FUN=mean)
SE<-aggregate(larvae~days*treatment*level, data=df.PyLv, FUN=std.error)

df.mean<-cbind(sum, SE[4])
colnames(df.mean)[4]<-"mean.larvae"
colnames(df.mean)[5]<-"SE"

larv.plot.mean<-ggplot(df.mean, aes(x=level, y=mean.larvae))+
  geom_point(aes(shape=treatment, color=treatment), position = position_dodge(width=0.9), size=3)+
  geom_errorbar(aes(color=treatment, ymin = mean.larvae - SE, ymax = mean.larvae + SE), 
                position = position_dodge(width=0.9), width=0, linewidth=0.8)+
  scale_color_manual(values=c("coral", "goldenrod", "plum4")) +
  facet_wrap(~treatment + days)+
  ylab("Surviving Larvae")+
  theme_bw()

larv.plot.mean
dev.copy(pdf, "figures/Larv.pyro.pdf", width = 6, height = 8)
dev.off()

df.PyLv.reg<-df.PyLv
df.PyLv.reg$days<-as.numeric(as.character(df.PyLv.reg$days))

larv.plot.reg<-ggplot(df.PyLv.reg, aes(x=days, y=larvae, fill=treatment))+
  facet_wrap(~treatment)+
  geom_jitter(aes(shape=level), size=5, alpha=0.5, width=0.1)+
  scale_fill_manual(values=c("dodgerblue3", "goldenrod", "coral")) +
  scale_color_manual(values=c("dodgerblue3", "goldenrod", "coral")) +
  scale_shape_manual(values=c(21,22,23,24)) +
  geom_smooth(method=lm, linewidth=2, aes(color=treatment))+
  ylab("Surviving Larvae")+
  theme_bw()

larv.plot.reg
dev.copy(pdf, "figures/Larv.pyro.regr.pdf", width = 10, height = 5)
dev.off()

larv.plot.mean
dev.copy(pdf, "figures/Larv.pyro.pdf", width = 6, height = 8)
dev.off()

mod.mix.lm<-lm(larvae~days+treatment, data=df.PyLv.reg)
check_model(mod.mix.lm)

mod.mix<-lmer(larvae~days*treatment + (1|level), data=df.PyLv.reg)
summary(mod.mix)
visreg::visreg(mod.mix, "days", by = "treatment")
Anova(mod.mix, type=2)
plot(allEffects(mod.mix))

mod.posthoc<-lmer(larvae~days+treatment + (1|level), data=df.PyLv.reg)
posthoc.mod<-emmeans(mod.posthoc, ~treatment)
multcomp::cld(posthoc.mod, Letters=LETTERS)

