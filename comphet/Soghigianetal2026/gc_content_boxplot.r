library(tidyverse)
library(cowplot)
library(ggtext)
library(viridis)
library(ggtext)

dat<-read.table('residues_per_tax_per_aln.txt', header = T, sep = '\t')
dat_out<-read.table('base_comp.txt', header = T)
dat_out<-dat_out %>%
  rename(tax = Sequence_ID, pos = Position, gene = Gene)

dat<-dat %>%
  rename(tax = Sequence_ID,
         gene_pos = File) %>%
  # head %>%
  mutate(gene_pos = gsub('(p[0-9]+_)(.+)', '\\2', gene_pos),
         tax = gsub(' ', '', tax)) %>%
  separate_wider_delim(cols = gene_pos, delim = '_', names = c('gene', 'pos'))

outgroup_first<-c('Mochlonyx', 'CLUMA', 'Chaoborus', 'Culicoides', 'Pvanderplanki') ##old set
old_tax<-unique(dat_out$tax)
old_out<-old_tax[grepl(paste(outgroup_first, collapse = "|"), old_tax)]
old_out<-old_out[!(old_out %in% outgroup)]

dat_out<-dat_out %>%
  select(gene, tax, pos, GC) %>%
  filter(tax %in% old_out) %>%
  mutate(pos = factor(pos, levels = c('1', '2', '3')))



outgroup_set<- c('Mochlonyx', 'Belgica', 'Chaoborus', 'Simulium', 'Forcipomyia') ##contains low and high GC outgroups
taxa<-levels(factor(dat$tax))
taxa<-gsub(' ', '', taxa)
outgroup<-taxa[grepl(paste(outgroup_set, collapse = "|"), taxa)]
ingroup<-taxa[!(taxa %in% outgroup)]
ingroup<-ingroup[ingroup != 'Wyeomyia_smithii']
ingroup<-ingroup[ingroup != 'Tamboinensis']
tox<-'Tamboinensis'
wyeo<-'Wyeomyia_smithii'
out_new<-c(old_out, outgroup)

dat_comb<-dat %>% 
  select(tax, gene, pos, C_Percent, G_Percent) %>%
  mutate(GC = C_Percent + G_Percent) %>%
  select(gene, tax, pos, GC) %>%
  bind_rows(dat_out) %>%
  mutate(grouping = case_when(tax %in% out_new ~ 'Outgroup',
                              tax %in% ingroup ~ 'Ingroup',
                              tax == 'Wyeomyia_smithii' ~ 'Wyeomyia smithii',
                              tax == tox ~ 'Toxorhynchites amboinensis'),
         labs = case_when(grouping == 'Wyeomyia smithii' ~ 'Sabethini 2<br>*Wyeomyia smithii*',
                          grouping =='Toxorhynchites amboinensis' ~ 'Toxorhynchitini<br>*Toxorhynchites amboinensis*',
                          grepl('Belgica', tax) ~ 'Outgroup 1<br>*Belgica antarctica*',
                          grepl('CLUMA', tax) ~ 'Outgroup 2<br>*Clunio marinus*',
                          grepl('Culicoides', tax) ~ 'Outgroup 3<br>*Culicoides sonorensis*',
                          grepl('Pvand', tax) ~ 'Outgroup 4<br>*Polypedilum vanderplanki*',
                          grepl('Moch', tax) ~ 'Outgroup 5<br>*Mochlonyx cinctipes*',
                          grepl('Chaob', tax) ~ 'Outgroup 6<br>*Chaoborus flavidulus*',
                          grepl('Forci', tax) ~ 'Outgroup 7<br>*Forcipomyia taiwana*',
                          grepl('Simuli', tax) ~ 'Outgroup 8<br>*Simulium* sp.',
                          .default = 'Ingroup'
         ),
         labs = factor(labs, levels = c('Outgroup 8<br>*Simulium* sp.',
                                        'Outgroup 7<br>*Forcipomyia taiwana*',
                                        'Outgroup 6<br>*Chaoborus flavidulus*',
                                        'Outgroup 5<br>*Mochlonyx cinctipes*',
                                        'Outgroup 4<br>*Polypedilum vanderplanki*',
                                        'Outgroup 3<br>*Culicoides sonorensis*',
                                        'Outgroup 2<br>*Clunio marinus*',
                                        'Outgroup 1<br>*Belgica antarctica*',
                                        'Toxorhynchitini<br>*Toxorhynchites amboinensis*',
                                        'Sabethini 2<br>*Wyeomyia smithii*',
                                        'Ingroup'))
  )


(p1<-dat_comb %>%
    # filter(grouping != 'Outgroup') %>%
    ggplot(aes(x = labs, y = GC, fill = as.factor(pos))) +
    geom_boxplot(position = position_dodge(width = 0.50),
                 width = 0.35, outliers = F, median.color = 'white') +
    # coord_cartesian(ylim = c(-100, 100)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(y = 'GC-content (%)', x = NULL) +
    coord_flip() +
    scale_fill_viridis(discrete = T, 'Codon position') +
    theme_bw() +
    theme(aspect.ratio = 1.75,
          legend.position = 'top',
          axis.text.y = element_markdown(color = 'black'),
          axis.text.x = element_text(color = 'black'),
          axis.title = element_text(color = 'black'))
)