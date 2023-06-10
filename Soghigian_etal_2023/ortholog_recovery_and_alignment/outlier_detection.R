
library(ape)
library(phytools)
library(plyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(adegenet)



#Set the working directory. This directory must have the following subdirs:
#a directory with trees to evaluate (in this example, trees_aa)
#a directory with nucleotide alignments (in this example, aln)
#a directory with amino acid alignments (in this example, aln_aa)

setwd("E:\\projects\\ahtree\\transcriptomes")

#create empty objects for the loop below.
per_tax_bl_ratio <- NULL
total_trees <- NULL

#loop over the input directory. Change the subdirectory and the extension as 
#needed to properly load in the trees for your particular dataset.
#Eventually, this script will be an R function that takes input of tree locations.

for (i in 1:length(list.files(path = (paste0(getwd(), "/trees_aa/")), pattern='*.treefile'))){
  loci <- list.files(path = (paste0(getwd(), "/trees_aa/")), pattern='*.treefile')[i]
  loci_name <- gsub(".treefile","",list.files(path=(paste0(getwd(), "/trees_aa/")), pattern='*.treefile')[i])
  tree <- read.tree((paste0(getwd(), "/trees_aa/",loci)))
  total_trees[[i]] <- tree
  bl_ratio <- as.data.frame(
    abs((colSums(cophenetic(tree))-median(colSums(cophenetic(tree))))/(IQR(colSums(cophenetic(tree)))+0.0000000000000000000001)
  ))
  colnames(bl_ratio)<-"bl_ratio"
  cbind(rownames(bl_ratio),loci_name,bl_ratio) -> per_tax_bl_ratio_tmp
  rbind(per_tax_bl_ratio,per_tax_bl_ratio_tmp) -> per_tax_bl_ratio
  rm(bl_ratio)
}

#Properly name our columns
colnames(per_tax_bl_ratio) <- c("tax","loci","bl_ratio")
#Let's make sure our trees are the right kind of object for plotting
class(total_trees) <- "multiPhylo"

#Define your outlier level. This is the number of interquartile ranges past the 
#IQR at which we consider something an outlier. Tukey said that 3 was extreme 
#outliers, but we'll be particularly conservativ in this case - and set it to 5.
#Note you can change this and evaluate how it performs below

outlierLevel <- 5
#What outliers exist? Let's evaluate our branch length ratios and find out
taxOutlier <- ddply( per_tax_bl_ratio, "tax", function(x) 
  data.frame( x, x$bl_ratio > (quantile(x$bl_ratio, 0.75) + (outlierLevel * IQR(x$bl_ratio))) ) )
colnames(taxOutlier)[4] <- 'isOutlier'
#Ok, how many did we find?
length(levels(as.factor(as.character(taxOutlier[taxOutlier$isOutlier==TRUE,]$loci))))
#Let's see how many genes have how many outliers
table(taxOutlier$loci,taxOutlier$isOutlier) -> loci_summaries

#We can print this object - or, we can save it for later:
write.csv(file="loci_summaries_5.csv",loci_summaries)

#We should also look at taxa summaries - these could be important in finding 
#problematic samples!
table(taxOutlier$tax,taxOutlier$isOutlier) -> tax_summaries

#We can print this object - or, we can save it for later:
write.csv(file="tax_summaries_5.csv",tax_summaries)



#We can trip alignments, assuming they are in subfolders /aln/ and /aln_aa/, as follows:

for (level in levels(as.factor(as.character(taxOutlier[taxOutlier$isOutlier==TRUE,]$loci)))){
  read.dna((paste0(getwd(), "/aln/",level,".fasta")), format="fasta") -> example1
  paste0(taxOutlier[with(taxOutlier,taxOutlier$loci==level & taxOutlier$isOutlier==TRUE),]$tax) -> to_remove
  ix <- which(rownames(example1) %in% to_remove)
  reduced_example1 <- example1[-ix, ]
  rm(example1)
  write.FASTA(reduced_example1, file=paste0(getwd(), "/aln/",level,"_pruned.fasta"))
  rm(reduced_example1)
  rm(to_remove)
  
  read.FASTA((paste0(getwd(), "/aln_aa/",level,".fasta")),type="AA")->aa_example1
  paste0(taxOutlier[with(taxOutlier,taxOutlier$loci==level & taxOutlier$isOutlier==TRUE),]$tax) -> to_remove
  ix <- which(names(aa_example1) %in% to_remove)
  reduced_aa_example1 <- aa_example1[-ix]
  rm(aa_example1)
  write.FASTA(reduced_aa_example1, file=paste0(getwd(), "/aln_aa/",level,"_pruned.fasta"))
  rm(reduced_aa_example1)
  rm(to_remove)
}

#We might want to seee the trimmed trees, potentially to plot untrimmed and trimmed trees, etc:
new_trees <- total_trees
#we need to set the loci list to factors, so they can be properly looped over.
as.factor(taxOutlier$loci) -> taxOutlier$loci
for (level in levels(as.factor(as.character(taxOutlier[taxOutlier$isOutlier==TRUE,]$loci)))){
  new_trees[[which(levels(taxOutlier$loci)==level)]] <- drop.tip(total_trees[[which(levels(taxOutlier$loci)==level)]],
                                                                 paste0(taxOutlier[with(taxOutlier,taxOutlier$loci==level & taxOutlier$isOutlier==TRUE),]$tax)) 
}


p1_hist1 <- ggplot(taxOutlier,aes(x=bl_ratio))+
  geom_histogram(bins=100,fill="goldenrod",color="black")+
  xlab("Per taxon per loci Branch Length Ratio")+
  ylab("Number of Loci")+
  theme_cowplot()

p2_box1 <- ggplot(taxOutlier,aes(x=tax,y=bl_ratio))+geom_boxplot(outlier.size=0.5)+
  xlab("Taxon")+ylab("Branch Length Ratio")+
  theme_cowplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=4))

taxOutlier[taxOutlier$isOutlier==FALSE,] -> notOutlier

p3_hist2 <- ggplot(notOutlier,aes(x=bl_ratio))+
  geom_histogram(bins=100,fill="goldenrod",color="black")+
  xlab("Post trim per taxon per loci Branch Length Ratio")+
  ylab("Number of Loci")+
  theme_cowplot()

p4_box2 <- ggplot(notOutlier,aes(x=tax,y=bl_ratio))+geom_boxplot(outlier.size=0.5)+
  xlab("Taxon")+ylab("Branch Length Ratio")+
  theme_cowplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=4))

four_panel <- plot_grid(p1_hist1,p2_box1,p3_hist2,p4_box2,labels = c('A', 'B', 'C','D'), label_size = 12)
save_plot(file="four_panel_gen_tra_t5.pdf",four_panel,base_height=8.5,base_width=11,scale=2)


#We can look at an individual locus... If we specify it here:
focus_locus="OG100"
focus_tree <- total_trees[[which(levels(as.factor(per_tax_bl_ratio$loci))==focus_locus)]]
outlierTax <- taxOutlier[with(taxOutlier,taxOutlier$loci==focus_locus & taxOutlier$isOutlier==TRUE),]$tax
tt<-paintBranches(focus_tree,edge=sapply(outlierTax,match,focus_tree$tip.label),
                  state="outlierTax",anc.state="a")
cols<-setNames(c("black","red"),c("a","outlierTax"))
plot(tt,colors=cols,lwd=2,split.vertical=TRUE,ftype="i")

#You can expect an individual species, too. In this case, we'll look at the example taxa pscil:
taxOutlier[taxOutlier$tax=="pscil",]
