library(tidyverse)
library(deeptime)
library(ggstance)
library(viridis)
library(cowplot)
dat<-read.table("hpd.txt", sep = "\t", header = TRUE)


tmp<-data.frame(start = c(0, 2.6, 23, 66, 143.1, 201.4, 251.9))
source_cols<-c('black', '#3F4889FF', '#1FA187FF', '#ed8911')

dat2<-dat %>%
  mutate(node = case_when(grepl('chaoboridae', node) ~ 'Chaoboridae', .default = node),
         node = factor(node, levels = c('Chaoboridae', 'Anophelinae', 'Culicidae', 'Culicinae',
                                        'Sabethini', 'Culisetini', 'Culicini', 'Aedini')),
         highlight = case_when(node %in% c('Aedini', 'Culisetini', 'Culicinae', 'Anophelinae') ~ TRUE, .default = FALSE),
         source = case_when(grepl('S2023', source) ~ 'SEA2023',
                            grepl('P2025', source) ~ 'PEA2025',
                            grepl('L2021', source) ~ 'LEA2021',
                            .default = source),
         fossil_tf = case_when(source == 'Fossil' ~ T, .default = F)) %>%
  mutate(taxa = fct_reorder(taxa, HPD_MID))


(p3<-ggplot(dat2, aes(x = HPD_MID, y = node, group = taxa, color = source)) +
    geom_tile(aes(x = 0, width = Inf, fill = highlight), height = 0.9, color = NA) +
    geom_segment(aes(x = HPD_MIN, xend = HPD_MAX, y = node, group = taxa, color = source),
                 position = position_dodgev(height = 0.75)) + #plot HPD interval, dodge each bar vertically
    geom_point(position = position_dodgev(height = 0.75), size = 2.5) + #plot HPD midpoint, dodge each vertically
    scale_x_reverse('Age (Ma)') + #set x-axis to go in reverse
    coord_geo(xlim = c(260, 0), ylim = c(0.60, 9)) + #set geoglogical time names
    scale_color_manual(values = source_cols, 'Source') +
    geom_vline(data = tmp, aes(xintercept = start), lwd = 0.25, color = 'grey78', lty = 2) +
    scale_fill_manual(values = c('grey88', 'white'), guide = 'none') +
    scale_shape_manual(values = c(16, 17), guide = 'none') +
    theme_classic() +
    theme(aspect.ratio = 12/8,
          axis.title.y = element_blank(),
          axis.title.x = element_text(color = 'black'),
          axis.text = element_text(color = 'black', size = 8)))