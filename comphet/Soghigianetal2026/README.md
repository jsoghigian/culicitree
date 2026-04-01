## Introduction
Here we provide the code used to generate the figures shown in Soghigian et al. 2026, "The Culicinae are Monophyletic and Ancient: A response to Pierce et al. 2025."  We provide a mix of commands used to generate files, as well as R scripts used to create visualizations. 

### Datasets  
- aa_ml_subgenus.tre - A phylogeny from Soghigian et al. 2023 
- hpd.txt - A tab delimited file of fossil ages and HPDs.
- base_comp.txt - Base composition ( GC content) values for Dataset A.
- residues_per_tax_aln.txt - Base composition (GC content) values for Dataset C.
- DatasetA_nt.phy - The original alignment of 54 genomes and transcriptomes published in Soghigian et al. 2023. Also has an associated nt.123.parts and nt.3.parts file, which are partition files.
- DatasetA_recodeRY.fasta and DatasetA_recode01.fasta  - The alignment of 54 genomes and transcriptomes published in Soghigian et al. 2023, but recoded following a 01/RY strategy. Also has an associated nt.123.parts and nt.3.parts file, which are partition files.
- DatasetC_net.fasta - The taxa from Dataset A but with 3 new outgroups: Forpomyia taiwanana - GCA_963930915.1 (Ceratapogonidae), Belgica antarctica - GCA_000775305.1 (Chironomidae), and a transcriptome from a Similium sp. - GGBP00000000 (Simuliidae). Also has an associated nt.123.parts and nt.3.parts file, which are partition files.

### Figure 1
![CulicidaePhylogeny](https://github.com/jsoghigian/culicitree/blob/main/comphet/images/fig1.png?raw=true)  
Figure 1 is based on previously published phylogenies and other data. The phylogeny, which we have uploaded here as aa_ml_subgenus.tre, is from [Soghigian et al. 2023](https://doi.org/10.1038/s41467-023-41764-y). This phylogeny was visualized in [FigTree](https://tree.bio.ed.ac.uk/software/figtree/), then modified in Illustrator to include unsampled lineages and to mark morphological characters. The other panel of this figure was generated in R. The scripts to do so are included in fossil_age_plot.R.


### Figure 2
![new OGS resolve NT3](https://github.com/jsoghigian/culicitree/blob/main/comphet/images/multi_panel_new_og.png?raw=true)  
We noticed substantial compositional heterogeneity in the form of GC content variation (Panel A), which we felt manifested as long branch attraction among ingroup species with low GC content and outgroup taxa with low GC content (Panel B) and a discrepancy between amino acids/nt1 and 2 and nucleotide position 3. Either recoding, or using different outgroups, resolved this discrepancy.  
#### GC Content
To visualize GC content, we used a boxplot.  We generated GC content values using [PhyKit](https://github.com/jlsteenwyk/phykit) for alignments with the following command:
```
cd aln/
for file in $(ls *.fasta);
do
og=$(basename $file .fasta)
phykit gc_content $file -v > gctmp.txt
paste gctmp.txt <(for line in $(awk -F'\t' '{print $2}' gctmp.txt);do echo ${og};done) >> ../gc.txt
rm gctmp.txt
done
```
The resulting GC values were analyzed in R, as described in the R script gc_content_boxplot.r . 
#### Phylogenetic estimation of NT3 topology
Dataset A was used to approximate Pierce et al. This dataset contains 54 genomes and transcriptomes and was previously published in Soghigian et al., but the NT3 position along was not analyzed. As such, we show the code for that analysis here:
```
iqtree3 -s DatasetA_nt.phy -spp DatasetA_nt.3.parts -pre nt3.MFP -m MFP -B 1000 -nt 16
```
An alternative would be to extract alignment columns, e.g. with AMAS, but because IQ-Tree is flexible one can exclude positions they do not wish to analyze, e.g. nucleotide positions 1 and 2, which is what we did here.
The resulting topology is shown above in panel A.
#### Recoding NT3
To address compositional heterogeneity at third codon position, we use an RY recoding scheme, wherein purines (A and G) are recoded to R and pyrimidines (T and C) are recoded to Y. This strategy can also be done in a binary fashion, in which R is 0 and Y is 1. Because neither RY nor 01 approaches yielded differences in topology, but because from a theoretical level analyzing recoded data as binary data has some justification (see IQTree help discussions on the topic), we present those codes here.  We used  [PhyKit](https://github.com/jlsteenwyk/phykit) for recoding:
```
phykit alignment_recoding DatasetA.fasta RY01-recoding.txt
```
This resulted in a recoded alignment, saved here as DatasetA_recode01. This dataset was then analyzed with IQTree:
```
iqtree3 -s DatasetA_recode01.fasta -spp DatasetA_recode01.3.parts -m GTR2+FO+R+I -nt 16 --prefix DatasetA_recodeRYas01
```
The resulting phylogeny is shown in panel B, above.

#### Analyzing New Outgroups  
We replaced the lowest GC content outgroups, *Clunio marinus* (Chironomidae), *Culicoides sonorensis* (Ceratapogonidae), and *Polypedilum vanderplanki* (Chironomidae), all of which have 3rd position GC content of 35% or less. Using the same methods described in Soghigian et  al. 2023, we used the protein sequences from *Forpomyia taiwanana* - GCA_963930915.1 (Ceratapogonidae), *Belgica antarctica* - GCA_000775305.1 (Chironomidae), and a transcriptome from a *Similium sp.* - GGBP00000000 (Simuliidae), and added these sequences to the alignments of the ingroup mosquitoes from Soghigian et al. 2023, resulting in alignments that differered in outgroup taxon only.  Again following methods in Soghigian et al. 2023, we then concatenated alignments and analyzed them together. 
For expediency, we specified the model GTR:
```
iqtree3 -s DatasetC_nt.fasta -spp DatasetC_nt.3.parts -pre s2026.base.nt3.gtrfo -m GTR+FO+R4 -B 1000 -nt 16
```
We later retrieved the same topology with ModelFinder.
