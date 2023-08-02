library(phytools)
library(dplyr)
library(tidyr)
library(tibble)
library(geiger)
library(geomorph)
library(ggplot2)
library(cowplot)
library(ggtree)
library(parallel)
library(plotrix)
library(foreach)
library(viridisLite)
library(doParallel)


#load host data set, taxonomy, and MCC of 100 TACT phylogenies
host.dat<-read.csv('species_blood_numerical_updatedJSS.csv', header=T, fill=T)
tax.dat<-read.csv('Culicidae_taxonomy_genusgroups.csv', header=T, fill=T)
mcc.tact.tree<-read.tree('mcc.100.tacted_071421.tre')

#format host data set for analysis
host_tax.dat<-tax.dat %>%
  select(genus.species, tribe) %>%
  mutate(genus.species = gsub(' ', '_', genus.species)) %>% #substitute all spaces with underscores ('_')
  right_join(., host.dat, by=c('genus.species' = 'Row.Labels')) %>% #match each species to their tribe
  column_to_rownames('genus.species') %>% #make the genus.species column into rownames
  select(tribe,amp_prop:rep_prop) %>% #subset the data set to include only columns necessary for analysis
  filter(tribe != 'Ficalbiini' & tribe != 'Hodgesiini') #Ficalbiini and Hodgesiini species have no molecular data, so remove them


tree.analysis<-treedata(mcc.tact.tree, host_tax.dat %>% select(!tribe))$phy #prune tree to match species in host_tax.dat
dat.analysis<-treedata(mcc.tact.tree, host_tax.dat %>% select(!tribe))$data #remove rows that are not in tree.analysis

#ancestral state reconstruction in parallel
#three models, equal rates (ER) symmetrical (SYM) and all rates different  (ARD)
#run each in parallel, 10 runs of a 100 simulations each.
#ER
t0<-Sys.time()
setCores<-round(detectCores()*0.85)
cl<-parallel::makeCluster(getOption('cl.cores', setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(phytools))
parallel::clusterExport(cl, 'tree.analysis')
parallel::clusterExport(cl, 'dat.analysis')
NRUNS <- 10
er.res.par<-parallel::parLapply(cl=cl, 1:NRUNS, function (x){
  make.simmap(tree.analysis, as.matrix(dat.analysis), model='ER', nsim=100)
})
stopCluster(cl)
t1<- Sys.time()
er.time<-t1-t0


#SYM
t0<-Sys.time()
setCores<-round(detectCores()*0.85)
cl<-parallel::makeCluster(getOption('cl.cores', setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(phytools))
parallel::clusterExport(cl, 'tree.analysis')
parallel::clusterExport(cl, 'dat.analysis')
NRUNS <- 10
sym.res.par<-parallel::parLapply(cl=cl, 1:NRUNS, function (x){
  make.simmap(tree.analysis, as.matrix(dat.analysis), model='SYM', nsim=100)
})
stopCluster(cl)
t1<- Sys.time()
sym.time<-t1-t0

#ARD
t0<-Sys.time()
setCores<-round(detectCores()*0.85)
cl<-parallel::makeCluster(getOption('cl.cores', setCores))
cl.pkg<-parallel::clusterEvalQ(cl, library(phytools))
parallel::clusterExport(cl, 'tree.analysis')
parallel::clusterExport(cl, 'dat.analysis')
NRUNS <- 10
ard.res.par<-parallel::parLapply(cl=cl, 1:NRUNS, function (x){
  make.simmap(tree.analysis, as.matrix(dat.analysis), model='ARD', nsim=100)
})
stopCluster(cl)
t1<- Sys.time()
ard.time<-t1-t0

#collect results
er.res.1000<-do.call(c, er.res.par)
sym.res.1000<-do.call(c, sym.res.par)
ard.res.1000<-do.call(c, ard.res.par)

#set class of the results
class(er.res.1000)<-c('multiSimmap', class(er.res.1000))
class(sym.res.1000)<-c('multiSimmap', class(sym.res.1000))
class(ard.res.1000)<-c('multiSimmap', class(ard.res.1000))

#summarize the results
er.res.1000.sum<-summary(er.res.1000, plot = F)
sym.res.1000.sum<-summary(sym.res.1000, plot = F)
ard.res.1000.sum<-summary(ard.res.1000, plot = F)

#compute AICc and weights compare ard, sym, and er models
aic<-function(logL,k) 2*k-2*logL
simmap.res<-list(er.res.1000, sym.res.1000, ard.res.1000)
names(simmap.res)<-c('ER', 'SYM', 'ARD')
simmap.aic.res<-mapply(aic,unlist(lapply(simmap.res, function(x) x[[1]]$logL)), c(1, 6, 12))
(aic.w(simmap.aic.res))


#plot results of each model
#set host colors
host<-setNames(viridis(4), c('amp_prop', 'ave_prop', 'mam_prop', 'rep_prop'))

#create dataframe of taxon names, tribe membership, and tribe label. Anopheles is not
#a tribe, but it should still be labelled.
dat.for.plot<-tax.dat %>%
  mutate(genus.species = gsub(' ', '_', genus.species)) %>%
  select(genus.species, tribe) %>%
  left_join(data.frame(dat.analysis) %>% rownames_to_column('genus.species') %>% select(genus.species), .) %>%
  mutate(tribe.label = case_when(tribe == 'None' ~ 'Anopheles',
                                 tribe != 'None' ~ tribe))

#time-scale axis and rings that represent different geologic eras
#to make a plot with rings as era cutoffs like Figure 2, we need to create a dummy
#plot
age.seq<-seq(0, 190, by = 10)
age.seq[which((age.seq/10)%%2>0)]<-NA
geo<-matrix(c(max(nodeHeights(tree.analysis)), 145, 145, 66, 66, 23, 23, 2.58, 2.58,0), byrow=T, ncol=2)
rownames(geo)<-c('Jurassic', 'Cretaceous','Paleogene', 'Neogene','Quaternary')
plot(tree.analysis)
circ.leg<-geo.legend(leg=geo, colors=c('grey78', 'white', 'grey78', 'white', 'grey78'))
r<-max(circ.leg$leg[,1]) -circ.leg$leg[,2]


#plot result of the ancestral state reconstructions
#ARD
#plot a foundation of the tree with results, but everything is kept transparent
#because all of the other elements will be plotted over it
# par(mai = c(2,2,2,2))
# pdf('sym.mcc.tact_071421.pdf', width=8 ,height=8)
plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5)
for(i in nrow(circ.leg$leg):1){
  color<-paste(strsplit(circ.leg$colors[i],'')[[1]][1:7], collapse='')
  draw.circle(0,0, radius=r[i], col=color, border='transparent')
}
plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5, add = T)
nodelabels(pie = sym.res.1000.sum$ace, piecol = host, cex = 0.35)
time.ax<-axis(1, pos = -6, at = seq(0, 190, by = 10), cex.axis = 1, labels = F)
text(rev(time.ax), rep(-15, length(time.ax)), labels = age.seq,  cex=0.6)
text(mean(time.ax), -22, 'Time (MYA)', cex=0.8)
arc.cladelabels(text = 'A', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedeomyiini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'E', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culisetini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'D', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Mansoniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'C', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Sabethini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'B', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Uranotaeniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Aedini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedini') %>% 
                                                   pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Anopheles', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Anopheles') %>% 
                                                      pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Culicini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culicini') %>% 
                                                     pull(genus.species)), mark.node = F)
text('A: Aedeomyiini',x=-200,y=166, cex=0.75, pos=4)
text('B: Uranotaeiini',x=-200,y=160, cex=0.75, pos=4)
text('C: Sabathini',x=-200,y=154, cex=0.75, pos=4)
text('D: Mansoniini',x=-200,y=148, cex=0.75, pos=4)
text('E: Culisetini',x=-200,y=142, cex=0.75, pos=4)
text('Model: SYM', x=-200, y=190, cex=1, pos=4)
legend(legend = names(host), x = -200, y = -140, cex = 0.75, bty = 'n', fill = host)
# dev.off()

# par(mai = c(2,2,2,2))
# pdf('er.mcc.tact_071421.pdf', width=8 ,height=8)
plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5)
for(i in nrow(circ.leg$leg):1){
  color<-paste(strsplit(circ.leg$colors[i],'')[[1]][1:7], collapse='')
  draw.circle(0,0, radius=r[i], col=color, border='transparent')
}
plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5, add = T)
nodelabels(pie = er.res.1000.sum$ace, piecol = host, cex = 0.35)
time.ax<-axis(1, pos = -6, at = seq(0, 190, by = 10), cex.axis = 1, labels = F)
text(rev(time.ax), rep(-15, length(time.ax)), labels = age.seq,  cex=0.6)
text(mean(time.ax), -22, 'Time (MYA)', cex=0.8)
arc.cladelabels(text = 'A', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedeomyiini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'E', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culisetini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'D', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Mansoniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'C', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Sabethini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'B', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Uranotaeniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Aedini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedini') %>% 
                                                   pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Anopheles', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Anopheles') %>% 
                                                      pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Culicini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culicini') %>% 
                                                     pull(genus.species)), mark.node = F)
text('A: Aedeomyiini',x=-200,y=166, cex=0.75, pos=4)
text('B: Uranotaeiini',x=-200,y=160, cex=0.75, pos=4)
text('C: Sabathini',x=-200,y=154, cex=0.75, pos=4)
text('D: Mansoniini',x=-200,y=148, cex=0.75, pos=4)
text('E: Culisetini',x=-200,y=142, cex=0.75, pos=4)
text('Model: ER', x=-200, y=190, cex=1, pos=4)
legend(legend = names(host), x = -200, y = -140, cex = 0.75, bty = 'n', fill = host)
# dev.off()


plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5)
for(i in nrow(circ.leg$leg):1){
  color<-paste(strsplit(circ.leg$colors[i],'')[[1]][1:7], collapse='')
  draw.circle(0,0, radius=r[i], col=color, border='transparent')
}
plotTree(tree.analysis, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5, add = T)
nodelabels(pie = ard.res.1000.sum$ace, piecol = host, cex = 0.35)
time.ax<-axis(1, pos = -6, at = seq(0, 190, by = 10), cex.axis = 1, labels = F)
text(rev(time.ax), rep(-15, length(time.ax)), labels = age.seq,  cex=0.6)
text(mean(time.ax), -22, 'Time (MYA)', cex=0.8)
arc.cladelabels(text = 'A', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedeomyiini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'E', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culisetini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'D', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Mansoniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'C', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Sabethini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'B', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Uranotaeniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Aedini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Aedini') %>% 
                                                   pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Anopheles', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Anopheles') %>% 
                                                      pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Culicini', node = findMRCA(tree.analysis, tips = filter(dat.for.plot, tribe.label == 'Culicini') %>% 
                                                     pull(genus.species)), mark.node = F)
text('A: Aedeomyiini',x=-200,y=166, cex=0.75, pos=4)
text('B: Uranotaeiini',x=-200,y=160, cex=0.75, pos=4)
text('C: Sabathini',x=-200,y=154, cex=0.75, pos=4)
text('D: Mansoniini',x=-200,y=148, cex=0.75, pos=4)
text('E: Culisetini',x=-200,y=142, cex=0.75, pos=4)
text('Model: ARD', x=-200, y=190, cex=1, pos=4)
legend(legend = names(host), x = -200, y = -140, cex = 0.75, bty = 'n', fill = host)
# dev.off()