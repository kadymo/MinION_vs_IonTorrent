#Add new variable for plotting - alpha diversity
variable1 = sample_data(data.1)$Platform
variable2 = sample_data(data.1)$Type
sample_data(data.1)$PT <-  mapply(paste0, variable1, sep = "_", variable2)

##Alpha diversity
std <- subset_samples(data.1, Type == "Standard")
samples <- subset_samples(data.1, Type == "Sample")

p <- plot_richness(samples, x="Pipeline" , color = "Platform", measures = c("Shannon"))
p3 <- plot_richness(std, x="Pipeline" , measures = c("Shannon"))+ facet_grid(~Database, scales = "free_x") + geom_jitter()
p2 = p + theme_bw() + geom_boxplot(data=p$data, aes(x=Pipeline, y=value, color=Platform), alpha=0.1) + 
  facet_grid(~Database, scales = "free_x") + geom_point(data=p$data, aes(x=Pipeline, y=value, color=Platform), position = position_dodge(width = 0.75)) +
  geom_point(data=p3$data, aes(x=Pipeline, y=value, shape=PT), colour = "black", position = position_dodge(width = 0.75))
p2$layers <- p2$layers[-1]
p2

#Alpha diversity testing
library(emmeans)
library(multcomp)
rich <- estimate_richness(samples, measures = c("Shannon"))
variable1 = sample_data(samples)$Platform
variable2 = sample_data(samples)$Pipeline
variable3 = sample_data(samples)$Database
rich$tab <-  mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3)
res.aov <- aov(Shannon ~ tab, data = rich)
summary(res.aov)
lsmeans <- lsmeans(res.aov, 'tab')

##should be adapted based on dataset
contrasts <- list(
  "ONT_r_pipe" =c(0,0,0,0,-1,0,1,0),
  "IT_r_pipe" =c(-1,0,1,0,0,0,0,0),
  "ONT_s_pipe" =c(0,0,0,0,0,-1,0,1),
  "IT_s_pipe" =c(0,-1,0,1,0,0,0,0),
  "Q_R_plat" =c(0,0,-1,0,0,0,1,0),
  "Q_S_plat" =c(0,0,0,-1,0,0,0,1),
  "M_R_plat" =c(-1,0,0,0,1,0,0,0),
  "M_S_plat" =c(0,-1,0,0,0,1,0,0),
  "ONT_Q_db" =c(0,0,0,0,0,0,-1,1),
  "IT_Q_db" =c(0,0,-1,1,0,0,0,0),
  "ONT_M_db" =c(0,0,0,0,-1,1,0,0),
  "IT_M_db" =c(-1,1,0,0,0,0,0,0)
)

sum_test <- summary(as.glht(contrast(lsmeans, contrasts)), test = adjusted("bonferroni"))
capture.output(sum_test, file = "sum_test.txt")

