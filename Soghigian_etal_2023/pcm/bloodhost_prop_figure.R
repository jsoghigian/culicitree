library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
library(tidyr)
library(ggdist)


dat$genus<-gsub('([A-z]+)_([a-z]+)','\\1', dat$Row.Labels)
dat$genus[296]<-'Culex'
dat$genus[70]<-'Aedes'
dat$genus[7]<-'Aedes'
dat$genus[226]<-'Culex'
dat$genus[363]<-'Limatus'
dat$genus[4]<-'Aedes'
dat$genus[108]<-'Anopheles'
dat$genus[313]<-'Culex'
dat$genus[316]<-'Culex'
dat$genus[377]<-'Mansonia'
dat$genus[421]<-'Uranotaenia'

tax_dat$genus.species<-gsub(' ', '_', tax_dat$genus.species)

dat2<-left_join(dat, tax_dat[c('genus.species', 'tribe')], by = c('Row.Labels' = 'genus.species')) %>%
  group_by(genus, tribe) %>%
  summarise(Amphibian = mean(amp_prop),
            Bird = mean(ave_prop),
            Mammal = mean(mam_prop),
            Reptile = mean(rep_prop))
dat2<-dat2[which(complete.cases(dat2)),]
dat2$tribe2<-NA
dat2$tribe2<-as.character(dat2$tribe2)
dat2[which(dat2$tribe == 'Aedeomyiini' | dat2$tribe == 'Aedini'), 'tribe2']<-'Aedeomyiini + Aedini'
dat2[which(dat2$tribe == 'Culicini' | dat2$tribe == 'Culisetini' | dat2$tribe == 'Ficalbiini' |
             dat2$tribe == 'Hodgesiini' | dat2$tribe == 'Uranotaeniini'), 'tribe2']<-'Culicini + Culisetini + Ficalbiini + Hodgesiini + Uranotaeniini'
dat2[which(dat2$tribe == 'Mansoniini' | dat2$tribe == 'None' | dat2$tribe == 'Sabethini'), 'tribe2']<-'Mansoniini + Anopheles + Sabethini'

dat2[which(dat2$tribe == 'None'), 'tribe']<-'Anophelinae'


dat2<- dat2 %>%
  pivot_longer(cols = c(Amphibian:Reptile), names_to = 'Host', values_to = 'props')

dat2[which(dat2$genus=='Uranotaenia'), 'props']<-dat2 %>% 
  filter(grepl('Urano', tribe)) %>%
  mutate(props = props * 1.065861) %>%
  pull(props)

dat2[which(dat2$genus=='Aedes'), 'props']<-dat2 %>% 
  filter(grepl('Aedes', genus)) %>%
  mutate(props = props * 1.010305)  %>%
  pull(props)

(p1<-ggplot(dat2, aes(fill = Host, y = props, x = tribe_genus)) +
  geom_bar(position = 'stack', stat = 'identity') +
  #facet_wrap(~ tribe, scales = 'free_x', nrow = 1) +
  theme_minimal_grid() +
  xlab('Genus') +
  ylab('Average Host Proportion') +
  scale_x_discrete(labels = dat2$genus, breaks = dat2$tribe_genus) +
  theme(aspect.ratio = 0.25,
        axis.text.x = element_text(angle = 90,
                                   face = 'italic',
                                   size = 10,
                                   vjust = 0.5,
                                   hjust = 0.95),
        strip.text = element_text(size = 11)) +
  scale_fill_viridis(discrete = T))

# ggsave2(filename = 'bloodhost_prop_bargraphs2.pdf', plot = p1,
# height = 7,
# width = 9,
# units = 'in')

host1<-setNames(viridis(4), c('Amphibian', 'Bird', 'Mammal', 'Reptile'))

host2<-setNames(viridis(4), c('amp_prop', 'ave_prop', 'mam_prop', 'rep_prop'))

(p2<-ggplot(dat2, aes(fill = Host, y = props, x = tribe_genus)) +
  geom_bar(position = 'stack', stat = 'identity') +
  facet_grid(. ~ tribe, scales = 'free_x', space = 'free') +
  theme_minimal_grid() +
  xlab('Genus') +
  ylab('Average Host Proportion') +
  scale_x_discrete(labels = dat2$genus, breaks = dat2$tribe_genus) +
  theme(axis.text.x = element_text(angle = 90,
                                   face = 'italic',
                                   size = 9,
                                   vjust = 0.5,
                                   hjust = 0.95),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 9)) +
  scale_fill_manual(values = host1))

# ggsave(file = 'bloodhost_prop_bargraphs4.pdf', plot = p2,
# width = 8,
# height = 5,
# units = 'in')

dat[296,1]<-'Culex_pipiensquinquefasciatus'

dat3<-dat %>%
  filter(Row.Labels != 'Culex_sp.' &
           Row.Labels != 'Aedes_' &
           Row.Labels != 'Culex_' &
           Row.Labels != 'Culex_spp' & 
           Row.Labels != 'Culex_spp.' &
           Row.Labels != 'Mansonia_sp.' &
           Row.Labels != 'Culiseta_spp' &
           Row.Labels != 'Ficalbia_spp' &
           Row.Labels != 'Anopheles_spp' &
           Row.Labels != 'Uranotaenia_sp.' & 
           Row.Labels != 'Aedes_baisasi' & 
           Row.Labels != 'Limatus_' &
           Row.Labels != 'Uranotaenia_sapphirina') %>%
  separate(Row.Labels, into = c('genus', 'species'), remove = F) %>%
  pivot_longer(cols = c(amp_prop:rep_prop), names_to = 'Host', values_to = 'prop') %>%
  group_by(Row.Labels) %>%
  summarise(prop = max(prop)) %>%
  left_join(., dat %>%
              filter(Row.Labels != 'Culex_sp.' &
                       Row.Labels != 'Aedes_' &
                       Row.Labels != 'Culex_' &
                       Row.Labels != 'Culex_spp' & 
                       Row.Labels != 'Culex_spp.' &
                       Row.Labels != 'Mansonia_sp.' &
                       Row.Labels != 'Culiseta_spp' &
                       Row.Labels != 'Ficalbia_spp' &
                       Row.Labels != 'Anopheles_spp' &
                       Row.Labels != 'Uranotaenia_sp.' & 
                       Row.Labels != 'Aedes_baisasi' & 
                       Row.Labels != 'Limatus_' &
                       Row.Labels != 'Uranotaenia_sapphirina') %>%
              pivot_longer(cols = c(amp_prop:rep_prop), names_to = 'Host', values_to = 'prop') %>%
              select(Row.Labels, genus, Host, prop)) %>%
  left_join(., tax_dat[c('genus.species', 'tribe')], by = c('Row.Labels' = 'genus.species'))

dat3[which(dat3$tribe == 'None'),'tribe']<-'Anophelinae'

library(ggdist)

(p3<-ggplot(dat3, aes(x = genus, y = prop, fill = Host)) +
  geom_jitter(shape =21, width = 0.15, alpha = 0.25) +
  facet_grid(.~tribe, space = 'free', scale = 'free_x') +
  theme_minimal_grid() +
  xlab('Genus') +
  ylab('Max Host Proportion') +
  scale_fill_manual(values = host2,
                    labels = c('Amphibians',
                               'Birds',
                               'Mammals',
                               'Reptiles')) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(angle = 90,
                                   face = 'italic',
                                   size = 9,
                                   vjust = 0.5,
                                   hjust = 0.95),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = 'bottom'))

# ggsave2(file = 'max_host_prop_all.pdf', plot = p3,
# width = 8,
# height = 5,
# units = 'in')

dat3_sub<-dat3 %>%
  filter(genus == 'Anopheles' |
           genus == 'Uranotaenia' |
           genus == 'Aedes' |
           genus == 'Psorophora' |
           genus == 'Culex' |
           genus == 'Culiseta' |
           genus == 'Coquillettidia' |
           genus == 'Mansonia')

(p4<-ggplot(dat3_sub, aes(x = genus, y = prop, fill = Host)) +
  geom_jitter(shape = 21, width = 0.15, alpha = 0.3) +
  theme_minimal_grid() +
  xlab('Genus') +
  ylab('Max Host Proportion') +
  scale_fill_manual(values = host2,
                    labels = c('Amphibians',
                               'Birds',
                               'Mammals',
                               'Reptiles')) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   face = 'italic',
                                   size = 9,
                                   vjust = 0.5,
                                   hjust = 0.95),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = 'bottom',
        aspect.ratio = 1.05))

# ggsave2(filename = 'max_host_prop_select_genera.pdf', plot = p4,
        # width = 8, height = 7, units = 'in')
