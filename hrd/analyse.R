#!/usr/bin/env Rscript
#
# Loads the file 'resultats-all.txt' as generated by the run-all.sh script and plots the values of all samples.
# Contains a list of samples known to be BRCA1/2 mutated or WT.
#
# The first argument is used to highlight a specific sample
		
args = commandArgs(trailingOnly=TRUE)
sample = args[1]

library(ggplot2)
library(gridExtra)

dat <- read.table("Q:/labo biol mol/Oncoscan/HRD-resultats_all.txt", sep='\t', header = TRUE)

vdat <- dat[,c('Nsegs','LOH','TAI','LST','TD','TDplus','Ploidy')]
vdat$diploid <- factor(ifelse(vdat$Ploidy<48,'diploid','non-diploid'))
vdat$mutations <- rep('N/A',dim(vdat)[1])
vdat$mutations[dat$DIAMIC=='H17003674'] <- 'BRCA1/2 mut�' #'BRCA2 VUS' #Baruchel, sein homme
vdat$mutations[dat$DIAMIC=='H17012299'] <- 'BRCA1/2 mut�' #'BRCA2 VUS' #Tavares, ovaire
vdat$mutations[dat$DIAMIC=='H18002311'] <- 'BRCA1/2 mut�' #'BRCA2 mut' #?
vdat$mutations[dat$DIAMIC=='H17015032'] <- 'BRCA1/2 mut�' #'BRCA1 VUS' #adca pulm
vdat$mutations[dat$DIAMIC=='H16013201'] <- 'BRCA1/2 mut�' #'BRCA1 mut' #ovaire
vdat$mutations[dat$DIAMIC=='H18000470'] <- 'BRCA1/2 mut�' #'BRCA1 mut' #Dery, ovaire
vdat$mutations[dat$DIAMIC=='H17005779'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17011883'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17013488'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17012993'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17012993-2'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17007513'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17013136'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17004324'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17013635'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17005622'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17007371'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H18002240'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H17012446'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H18005705'] <- 'BRCA1/2 WT'
vdat$mutations[dat$DIAMIC=='H18002163'] <- 'BRCA1/2 WT'

vdat$mutations[dat$DIAMIC=='H17010125'] <- 'CDK12 mut'
vdat$mutations[dat$DIAMIC=='H17000903'] <- 'BRCA1/2 mut�'
vdat$mutations[dat$DIAMIC=='H14002865'] <- 'BRCA1/2 mut�'
vdat$mutations[dat$DIAMIC=='H17002330'] <- 'BRCA1/2 mut�'
vdat$mutations[dat$DIAMIC=='H15012937'] <- 'BRCA1/2 mut�'
vdat$mutations[dat$DIAMIC=='H12003857'] <- 'BRCA1/2 mut�' #vrai DIAMIC: H12003845
vdat$mutations[dat$DIAMIC=='H11000871'] <- 'BRCA1/2 mut�'

vdat$mutations <- factor(vdat$mutations, 
                         levels = c('N/A','BRCA1/2 mut�','BRCA1/2 WT','CDK12 mut'))

vdat$selection <- rep('',dim(vdat)[1])
vdat$selection[dat$DIAMIC==sample] <- 'Sample'

#Make graphs
g.nsegs <- ggplot(data=vdat, aes(x=Ploidy, y=Nsegs)) +
  geom_jitter(width = 0.2, aes(col=mutations, shape=selection)) + ggtitle('Nombre de segments')
g.ploVSloh <- ggplot(data=vdat, aes(x=Ploidy, y=LOH)) + 
  geom_point(aes(col=mutations, shape=selection)) + ggtitle('Perte d\'h�t�rozygotie')
g.ploVStai <- ggplot(data=vdat, aes(x=Ploidy, y=TAI)) + 
  geom_point(aes(col=mutations, shape=selection)) + ggtitle('Alt�rations t�lom�riques')

g.ploVSlst <- ggplot(data=vdat, aes(x=Ploidy, y=LST)) + 
  geom_point(aes(col=mutations, shape=selection)) + ggtitle('Large state transitions')
g.ploVStdp <- ggplot(data=vdat, aes(x=Ploidy, y=TDplus)) + 
  geom_point(aes(col=mutations, shape=selection)) + ggtitle('Large tandem repeats')
g.ploVStd <- ggplot(data=vdat, aes(x=Ploidy, y=TD)) + 
  geom_point(aes(col=mutations, shape=selection)) + ggtitle('Small tandem repeats')
grid.arrange(g.nsegs, g.ploVSloh, g.ploVStai, g.ploVSlst, g.ploVStdp, g.ploVStd, ncol=3)


png(paste0('Q:/Biologie moleculaire/TB Molec/presentations/profilCNV-',sample,'.png'), width = 900, height = 700)
grid.arrange(g.nsegs, g.ploVSloh, g.ploVStai, g.ploVSlst, g.ploVStdp, g.ploVStd, ncol=3)
dev.off()



