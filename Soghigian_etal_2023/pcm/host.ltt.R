library(phytools)
library(dplyr)
library(tidyr)
library(tibble)
library(geiger)
library(ggplot2)
library(viridis)
library(deeptime)
library(cowplot)
library(ggnewscale)

host<-setNames(viridis(4), c('amp_prop', 'ave_prop', 'mam_prop', 'rep_prop'))

amph.tax<-read.csv('amph_shl_new_Classification.csv', header = T)
amph.tree<-read.nexus('amphibian_posterior.nex')
rep.tax<-read.csv('squamate_names.csv', header = T)
rep.tree<-read.nexus('squamata_posterior.nex')
aves.tree<-read.nexus('aves_posterior.nex')
aves.tax<-read.csv('birds_BLIOCPhyloMasterTax.csv', header = T)
mamm.tree<-read.nexus('mammalia_posterior.nex')
mamm.tax<-read.csv('taxonomy_mamPhy_5911species.csv', header = T)

amph.tax<-amph.tax %>%
  group_by(Family) %>%
  slice_head(n = 1) %>%
  mutate(Scientific.Name = gsub(' ', '_', Scientific.Name)) %>%
  column_to_rownames('Scientific.Name')

rep.tax<-rep.tax %>%
  group_by(Family) %>%
  slice_head(n = 1) %>%
  mutate(Species = gsub(' ', '_', Species)) %>%
  column_to_rownames('Species')

mamm.tax<-mamm.tax %>% 
  group_by(fam) %>%
  slice_head(n = 1) %>%
  column_to_rownames('Species_Name')

aves.tax<-aves.tax %>% 
  group_by(BLFamilyLatin) %>%
  slice_head(n = 1) %>%
  column_to_rownames('TipLabel')

class(amph.tree)<-c('multiPhylo', 'list')
amph.fam.trees<-lapply(amph.tree, function(x) treedata(x, amph.tax, sort = T, warnings = F)$phy)
class(amph.fam.trees)<-'multiPhylo'
amph.fam.ltt<-ltt95(amph.fam.trees, log = F, method = 'lineages', mode = 'median')
amph.fam.ltt.dat<-data.frame(amph.fam.ltt[, 1:4])
amph.fam.ltt.dat$new.lineages<-c(0, diff(amph.fam.ltt.dat$lineages))

class(rep.tree)<-c('multiPhylo', 'list')
rep.fam.trees<-lapply(rep.tree, function(x) treedata(x, rep.tax, sort = T, warnings = F)$phy)
class(rep.fam.trees)<-'multiPhylo'
rep.fam.ltt<-ltt95(rep.fam.trees, log = F, method = 'lineages', mode = 'median')
rep.fam.ltt.dat<-data.frame(rep.fam.ltt[, 1:4])
rep.fam.ltt.dat$new.lineages<-c(0, diff(rep.fam.ltt.dat$lineages))

class(mamm.tree)<-c('multiPhylo', 'list')
mamm.fam.trees<-lapply(mamm.tree, function(x) treedata(x, mamm.tax, sort = T, warnings = F)$phy)
class(mamm.fam.trees)<-'multiPhylo'
mamm.fam.ltt<-ltt95(mamm.fam.trees, log = F, method = 'lineages', mode = 'median')
mamm.fam.ltt.dat<-data.frame(mamm.fam.ltt[, 1:4])
mamm.fam.ltt.dat$new.lineages<-c(0, diff(mamm.fam.ltt.dat$lineages))

class(aves.tree)<-c('multiPhylo', 'list')
aves.fam.trees<-lapply(aves.tree, function(x) treedata(x, aves.tax, sort = T, warnings = F)$phy)
class(aves.fam.trees)<-'multiPhylo'
aves.fam.ltt<-ltt95(aves.fam.trees, log = F, method = 'lineages', mode = 'median')
aves.fam.ltt.dat<-data.frame(aves.fam.ltt[, 1:4])
aves.fam.ltt.dat$new.lineages<-c(0, diff(aves.fam.ltt.dat$lineages))


write.tree(amph.fam.trees, 'amph.fam.posterior.tre')
write.tree(aves.fam.trees, 'aves.fam.posterior.tre')
write.tree(mamm.fam.trees, 'mamm.fam.posterior.tre')
write.tree(rep.fam.trees, 'rep.fam.posterior.tre')

(all.classes.lines2<-ggplot() +
  stat_smooth(data = mamm.fam.ltt.dat[-101,], aes(x = abs(time - max(time)), y =  new.lineages), size = 1.5, color = '#35B779FF', se = F, span = 0.25) +
  new_scale_color() +
  stat_smooth(data = aves.fam.ltt.dat[-101,], aes(x = abs(time - max(time)), y = new.lineages), size = 1.5, color = '#31688EFF', se = F, span = 0.25) +
  new_scale_color() +
  stat_smooth(data = rep.fam.ltt.dat[-101,], aes(x = abs(time - max(time)), y =  new.lineages), size = 1.5, color = '#FDE725FF', se = F, span = 0.25) +
  new_scale_color() +
  stat_smooth(data = amph.fam.ltt.dat[-101,], aes(x = abs(time - max(time)), y =  new.lineages), size = 1.5, color = '#440154FF', se = F, span = 0.25) +
  coord_geo(xlim = c(180, 0)) +
  scale_x_reverse('Time (MYA)') +
  ylab('Number of new families') +
  theme_bw() +
  theme(aspect.ratio = 0.25))
# ggsave2(filename = 'all.classes.lines2.pdf', plot = all.classes.lines2, 
          # height = 4, width = 8, units = 'in')
