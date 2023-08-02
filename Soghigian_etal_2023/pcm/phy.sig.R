library(phytools)
library(dplyr)
library(tidyr)
library(tibble)
library(geiger)
library(geomorph)
library(phylosignal)
library(phylobase)

tree<-read.tree('mcc.100.tacted_071421.tre')
tree<-force.ultrametric(tree)
dat<-read.csv('species_blood_numerical_updatedJSS.csv', header=T, fill=T)
colnames(dat)[1]<-'genus.species'
tax.dat<-read.csv('Culicidae_taxonomy_genusgroups.csv', header=T, fill=T)
tax.dat$genus.species<-gsub(' ', '_',tax.dat$genus.species)



dat.analysis<-dat %>%
  left_join(tax.dat %>% select(genus.species, tribe)) %>%
  select(genus.species, tribe, amp_prop:rep_prop) %>%
  filter(tribe != 'Ficalbiini' &
           tribe != 'Hodgesiini') %>%
  select(!tribe) %>%
  column_to_rownames('genus.species')

dat.analysis<-data.frametreedata(phy = tree, data = dat.analysis)$data
tree.analysis<-treedata(phy = tree, data = dat.analysis)$phy


#estimate phylogenetic signal (K-mult) for all host classes together as well as each
#host class individually (K)
kmult<-physignal(as.matrix(dat.analysis), phy = tree.analysis, iter = 10000)
kmam<-physignal(as.matrix(select(dat.analysis, mam_prop)), phy = tree.analysis, iter = 10000)
kaves<-physignal(as.matrix(select(dat.analysis, ave_prop)), phy = tree.analysis, iter = 10000)
kamp<-physignal(as.matrix(select(dat.analysis, amp_prop)), phy = tree.analysis, iter = 10000)
krep<-physignal(as.matrix(select(dat.analysis, rep_prop)), phy = tree.analysis, iter = 10000)

#use phylosignal to estimate phylogenetic signal using a phylogenetic correlogram
#these take a while to run...
phylosignal.dat<-phylo4d(tree.analysis, dat.analysis)
amp_corr<-phyloCorrelogram(p4d = phylosignal.dat, trait = 'amp_prop')
ave_corr<-phyloCorrelogram(p4d = phylosignal.dat, trait = 'ave_prop')
rep_corr<-phyloCorrelogram(p4d = phylosignal.dat, trait = 'rep_prop')
mam_corr<-phyloCorrelogram(p4d = phylosignal.dat, trait = 'mam_prop')
all_corr<-phyloCorrelogram(p4d = phylosignal.dat, trait = c('amp_prop', 'ave_prop',
                                                            'rep_prop', 'mam_prop'))

#repeat above analyses without anophelines
dat.analysis2<-dat %>%
  left_join(tax.dat %>% select(genus.species, tribe)) %>%
  select(genus.species, tribe, amp_prop:rep_prop) %>%
  filter(tribe != 'Ficalbiini' &
           tribe != 'Hodgesiini' &
           tribe != 'None') %>%
  select(!tribe) %>%
  column_to_rownames('genus.species')

tree.analysis2<-treedata(tree.analysis, dat.analysis2)$phy

#estimate phylogenetic signal (K-mult) for all host classes as well as each
#host class individually (K)
kmult.res2<-physignal(as.matrix(dat.analysis2), phy = tree.analysis2, iter = 10000)
kmam2<-physignal(as.matrix(select(dat.analysis2, mam_prop)), phy = tree.analysis2, iter = 10000)
kaves2<-physignal(as.matrix(select(dat.analysis2, ave_prop)), phy = tree.analysis2, iter = 10000)
krep2<-physignal(as.matrix(select(dat.analysis2, rep_prop)), phy = tree.analysis2, iter = 10000)
kamp2<-physignal(as.matrix(select(dat.analysis2, amp_prop)), phy = tree.analysis2, iter = 10000)

bind_rows(mapply(function(x, y) {
  data.frame(host = y, phy.sig = x$phy.signal, x$pvalue, x$Z)}, 
  x = list(kmult.res2, kmam2, kaves2, krep2, kamp2), 
  y = list('all', 'mam', 'aves', 'rep', 'amp'), SIMPLIFY = F))

#estimate phylogenetic signal using a phylogenetic correlogram
phylosignal.dat2<-phylo4d(tree.analysis2, dat.analysis2)
amp_corr2<-phyloCorrelogram(p4d = phylosignal.dat2, trait = 'amp_prop')
ave_corr2<-phyloCorrelogram(p4d = phylosignal.dat2, trait = 'ave_prop')
rep_corr2<-phyloCorrelogram(p4d = phylosignal.dat2, trait = 'rep_prop')
mam_corr2<-phyloCorrelogram(p4d = phylosignal.dat2, trait = 'mam_prop')
all_corr2<-phyloCorrelogram(p4d = phylosignal.dat2, trait = c('amp_prop', 'ave_prop',
                                                            'rep_prop', 'mam_prop'))

# par(pty = 's', mai = c(2,2,2,2))
# pdf('amp_correlogram2.pdf', width = 6, height = 6)
plot(amp_corr2, main = 'Moran`s I for amphibian host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('ave_correlogram2.pdf', width = 6, height = 6)
plot(ave_corr2, main = 'Moran`s I for bird host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('mam_correlogram2.pdf', width = 6, height = 6)
plot(mam_corr2, main = 'Moran`s I for mammal host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('rep_correlogram2.pdf', width = 6, height = 6)
plot(rep_corr2, main = 'Moran`s I for reptile host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('all_correlogram2.pdf', width = 6, height = 6)
plot(rep_corr2, main = 'Mantel R for all hosts', cex.main = 0.75)
# dev.off()


# par(pty = 's', mai = c(2,2,2,2))
# pdf('amp_correlogram.pdf', width = 6, height = 6)
plot(amp_corr2, main = 'Moran`s I for amphibian host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('ave_correlogram.pdf', width = 6, height = 6)
plot(ave_corr2, main = 'Moran`s I for bird host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('mam_correlogram.pdf', width = 6, height = 6)
plot(mam_corr2, main = 'Moran`s I for mammal host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('rep_correlogram.pdf', width = 6, height = 6)
plot(rep_corr2, main = 'Moran`s I for reptile host %', cex.main = 0.75)
# dev.off()

# par(pty = 's', mai = c(2,2,2,2))
# pdf('all_correlogram.pdf', width = 6, height = 6)
plot(rep_corr2, main = 'Mantel R for all hosts', cex.main = 0.75)
# dev.off()