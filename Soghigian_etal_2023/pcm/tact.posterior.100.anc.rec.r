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
library(doParallel)
library(viridisLite)

host.dat<-read.csv('species_blood_numerical_updatedJSS.csv', header=T, fill=T)
tax.dat<-read.csv('Culicidae_taxonomy_genusgroups.csv', header=T, fill=T)
tact.tree<-read.tree('culicidae.tacted.100_071421.trees')
mcc.tact.tree<-read.tree('mcc.100.tacted_071421.tre')

dat.analysis<-tax.dat %>%
	select(genus.species, tribe) %>%
	mutate(genus.species = gsub(' ', '_', genus.species)) %>%
	right_join(host.dat, by = c('genus.species' = 'Row.Labels')) %>%
	# filter(tribe != 'Ficalbiini' &
		   # tribe != 'Hodgesiini') %>%
	column_to_rownames('genus.species') %>%
	select(amp_prop:rep_prop)


tact.tree.analysis<-lapply(tact.tree, function (x) treedata(x, dat.analysis)$phy)
dat.analysis<-treedata(tact.tree.analysis[[1]], dat.analysis)$data

t0<-Sys.time()
setCores<-round(detectCores()*0.85)
cl<-parallel::makeCluster(getOption('cl.cores', setCores))
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, 'tact.tree.analysis')
parallel::clusterExport(cl, 'dat.analysis')
tact.fit<-foreach(i=1:100, .packages = 'phytools') %dopar% {
  fitMk(tact.tree.analysis[[i]], dat.analysis, model = 'ARD')
}
stopCluster(cl)
t1<-Sys.time()  
t1-t0

fitQ.list<-lapply(tact.fit, function(x) {
  fitQ<-matrix(NA, length(x$states), length (x$states))
  fitQ[]<-c(0, x$rates)[x$index.matrix + 1]
  diag(fitQ)<-0
  diag(fitQ)<--rowSums(fitQ)
  colnames(fitQ)<-rownames(fitQ)<-x$states
  fitQ
})

t0<-Sys.time()
setCores<-round(detectCores()*0.85)
cl<-parallel::makeCluster(getOption('cl.cores', setCores))
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, 'tact.tree.analysis')
parallel::clusterExport(cl, 'dat.analysis')
parallel::clusterExport(cl, 'fitQ.list')
tact.ard.res<-foreach(i=1:100, .packages = 'phytools') %dopar% {
  make.simmap(tact.tree.analysis[[i]], dat.analysis, Q = fitQ.list[[i]], nsim = 100)
}
stopCluster(cl)
t1<-Sys.time()
t1-t0

tact.ard.each.summary<-lapply(tact.ard.res, summary, plot = F)
bind_rows(lapply(tact.ard.each.summary, function (x) data.frame(prob=x$ace[1,]) %>%
                   rownames_to_column('host') %>%
                   pivot_wider(names_from = host, values_from = prob)))
				   
mcc.tact.ref.tree<-treedata(mcc.tact.tree, dat.analysis)$phy

tact.ard.res.10000<-do.call(c, tact.ard.res)
class(tact.ard.res.10000)<- c('multiSimmap', class(tact.ard.res.10000))
t0<-Sys.time()
tact.ard.res.10000.summary<-summary(tact.ard.res.10000,  plot=F)
t1<-Sys.time()
t1-t0


host<-setNames(viridis(4), c('amp_prop', 'ave_prop', 'mam_prop', 'rep_prop'))

dat.for.plot<-tax.dat %>%
  mutate(genus.species = gsub(' ', '_', genus.species)) %>%
  select(genus.species, tribe) %>%
  left_join(data.frame(dat.analysis) %>% rownames_to_column('genus.species') %>% select(genus.species), .) %>%
  mutate(tribe.label = case_when(tribe == 'None' ~ 'Anopheles',
                                 tribe != 'None' ~ tribe))
mcc.tact.ref.tree<-rotateNodes(mcc.tact.ref.tree, findMRCA(mcc.tact.ref.tree, tips = c('Hodgesia_cyptopus', 'Ficalbia_ichiromiyagii')))

age.seq<-seq(0, 190, by = 10)
age.seq[which((age.seq/10)%%2>0)]<-NA
geo<-matrix(c(max(nodeHeights(mcc.tact.ref.tree)), 145, 145, 66, 66, 23, 23, 2.58, 2.58,0), byrow=T, ncol=2)
rownames(geo)<-c('Jurassic', 'Cretaceous','Paleogene', 'Neogene','Quaternary')
plot(mcc.tact.ref.tree)
circ.leg<-geo.legend(leg=geo, colors=c('grey78', 'white', 'grey78', 'white', 'grey78'))
r<-max(circ.leg$leg[,1]) -circ.leg$leg[,2]

# par(mai = c(2,2,2,2))
# pdf('ard.tact.posterior.10k_071421.pdf', width=8 ,height=8)
plotTree(mcc.tact.ref.tree, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5)
for(i in nrow(circ.leg$leg):1){
  color<-paste(strsplit(circ.leg$colors[i],'')[[1]][1:7], collapse='')
  draw.circle(0,0, radius=r[i], col=color, border='transparent')
}
plotTree(mcc.tact.ref.tree, ftype = 'off', type = 'fan', part = 0.97, lwd = 1.5, add = T)
nodelabels(pie = tact.ard.res.10000.summary$ace, piecol = host, cex = 0.35)
time.ax<-axis(1, pos = -6, at = seq(0, 190, by = 10), cex.axis = 1, labels = F)
text(rev(time.ax), rep(-15, length(time.ax)), labels = age.seq,  cex=0.6)
text(mean(time.ax), -22, 'Time (MYA)', cex=0.8)
arc.cladelabels(text = 'A', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Aedeomyiini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'F', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Culisetini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'E', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Mansoniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'C', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Sabethini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'B', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Uranotaeniini') %>% 
                                              pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Aedini', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Aedini') %>% 
                                                   pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Anopheles', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Anopheles') %>% 
                                                      pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'Culicini', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Culicini') %>% 
                                                     pull(genus.species)), mark.node = F)
arc.cladelabels(text = 'D', node = findMRCA(mcc.tact.ref.tree, tips = filter(dat.for.plot, tribe.label == 'Ficalbiini') %>% 
                                                     pull(genus.species)), mark.node = F)
add.arrow(tree = mcc.tact.ref.tree, tip = 'Hodgesia_cyptopus', col = 'black', lwd =1.5, arrl = 15)
text('A: Aedeomyiini',x=-200,y=166, cex=0.75, pos=4)
text('B: Uranotaeiini',x=-200,y=160, cex=0.75, pos=4)
text('C: Sabathini',x=-200,y=154, cex=0.75, pos=4)
text('D: Ficalbiini',x=-200,y=148, cex=0.75, pos=4)
text('E: Mansoniini',x=-200,y=142, cex=0.75, pos=4)
text('F: Culisetini',x=-200,y=136, cex=0.75, pos=4)
text('Model: ARD', x=-200, y=190, cex=1, pos=4)
legend(legend = names(host), x = -200, y = -140, cex = 0.75, bty = 'n', fill = host)
# dev.off()



plot(tact.ard.res.10000.summary2, fsize = 0.35, type ='fan', part = 0.97, lwd = 1.5, colors = host, cex = c(0.35,0.25))
text('Majority Rule-no reference specified', x = -200, y= 220)
