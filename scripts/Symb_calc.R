### symbiont fire-water larvae

symbs<-read.csv("data/cell counts.csv")
metad<-read.csv("data/firewater larvae_metadata.csv")

# merge dataframes
symb.df<- symbs %>%
  group_by(Sample.ID) %>%
  left_join(metadata, by="Sample.ID")

###### set structure
symb.df <- symb.df %>% 
  rename_with(~ str_remove(.x, "\\.x$"))

make.fac <- c("Sample.ID", "Treatment", "Level", "Rep")
symb.df[make.fac] <- lapply(symb.df[make.fac], factor)
symb.df$Treatment<- factor(symb.df$Treatment, levels=c("control", "burned", "inorganic"))
symb.df$Level<- factor(symb.df$Level, levels=c("control", "low", "medium", "high"))


# plot data
symb.dens.plot<- ggplot(symb.df, aes(interaction(Treatment, Level, sep = "\n "), cells.larvae, fill=Treatment))+
  geom_jitter(width=0.1, aes(color=Treatment), alpha=0.7)+
  geom_boxplot(outlier.shape = NA, alpha=0.6)+
  scale_fill_manual(values=c("dodgerblue", "coral", "goldenrod")) +
  scale_color_manual(values=c("dodgerblue", "coral", "goldenrod")) +
  ylab(expression(paste("symbiont density(cells larva"^{-1},")")))+
  scale_y_continuous(limits= c(1500, 4500))+
  xlab("Treatment State")+ theme_bw() + guides(color="none")

symb.dens.plot
dev.print(pdf, "figures/pyrolarv_symb.pdf", width = 7, height = 5)
dev.off()
