library(ggplot2, gridExtra)

sample <- c(208, 142, 55, 85)
dat.random <- read.table('Q:/labo biol mol/Oncoscan/results-100random_comparisons.txt', header = TRUE, sep = '\t')
dat <- rbind(dat.random, sample)

dat.norm <- data.frame(apply(dat, 2, function(v) v/rowSums(dat)))

dat.norm$type <- c(rep('random', dim(dat)[1]-1), 'sample')

ggplot(data=dat.melt, aes(x=variable, y=value)) + geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width= 0.2, aes(col=type))
qplot(log2(dat.norm$AtoB/dat.norm$BtoA), log2(dat.norm$AeqB/dat.norm$AneB), col=dat.norm$type)