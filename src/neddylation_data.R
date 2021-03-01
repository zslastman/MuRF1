# # Nedylation proteomic data comparison
#
# Background
# Data are collected from: https://www.ebi.ac.uk/pride/archive/projects/PXD022812
# The paper is [here](https://www.cell.com/cell-reports/fulltext/S2211-1247(20)31624-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124720316247%3Fshowall%3Dtrue).
#
# NEDD8 resembles ubiquitin. It is added by neddylases (canonical) and the ubiquitin pathway, (non canonical).
#
# The authors collected proteoimcs from 2 samples with denedylases knocked out, and 2 samples under proteoxic stress.
#

# ## load our existing fold changes
#
# We want to:
# - load the neddylation data
# - process it so it matches our fold changes

library(tidyverse);
library(data.table);
# library(ggExtra)
library(magrittr);

# +
library(here)
#readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3)
mg132_gnms = suppressMessages(readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3))[[2]]%>%.[!is.na(.)]
nedd8ko_gnms = suppressMessages(readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3))[[5]]

length(mg132_gnms)
length(nedd8ko_gnms)

mg132_gnms%>%tail
nedd8ko_gnms%>%tail

lobato_gnms = bind_rows(.id='lobato_set', mg132=tibble(gnm=mg132_gnms), nedd8ko = tibble(gnm=nedd8ko_gnms))
tail(lobato_gnms)
#readxl::excel_sheets('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3)
# -

# So it looks like we have 4 samples here, the MG_1 and 2, and NEDP1 and 2. I'll assume these are the proteotoxic and NEDD8 KO datasets respectively. Since there's no control group to compare to (it seems), we can just averag eintensities in our 

# ## Harmonize gene identifiers
#
# We'll use a file downloaded from the biomart website to map our human gene names -> human ids -> mouse ids -> mouse names
#

# +
orthodf = fread(here('ext_data/human_mouse_orthologs_biomart.txt.gz'))
hgnmdf = fread('../ext_data/human_gid_gnm_biomart.txt.gz')

hgnmdf = hgnmdf%>%gather(src,gnm,`Gene Synonym`:`Gene name`)%>%distinct(`Gene stable ID`,`gnm`)

#hgnmdf%<>%rename('gnm':=`Gene name`)
# hgnmdf%>%head

message('genes for which we can get the ID from lobatos gene names')
message('ned8 genes')
table(nedd8ko_gnms %in% hgnmdf$gnm)
message('mg132 genes')
table(mg132_gnms %in% hgnmdf$gnm)

lobato_gnm_orth = lobato_gnms%>%inner_join(hgnmdf,'gnm')%>%left_join(orthodf,by="Gene stable ID")
# lobato_gnm_orth%>%head
lobato_gnm_orth%<>%rename('mgid':=`Mouse gene stable ID`)


message('for how many of our combined lobato gene names do we  find an ortholog?')
lobato_gnm_orth%>%group_by(gnm)%>%summarise(noorth=any(mgid!=''))%>%.$noorth%>%table

# -

# # Now we load up our fold changes as well

resultslist <- readRDS(here('data/resultslist.rds'))
mrna_gids = resultslist[[1]]$gene_id
message('ow many of these genes have at least one ortholog appearing our annotation?')
lobato_gnm_orth%>%group_by(gnm)%>%summarise(inmrna_dat = any(mgid %in% mrna_gids))%>%.$inmrna_dat%>%table


# +
options(jupyter.plot_mimetypes = "image/svg+xml") 



contrasts = names(resultslist)
contrast = contrasts[[1]]

lobato_lfc_df = resultslist[[contrast]]%>%left_join(lobato_gnm_orth,by=c('gene_id'='mgid'))%>%distinct(lobato_set,gnm,log2FoldChange)
lobato_lfc_df%>%sample_n(4)

lobato_lfc_df%>%qplot(data=.,fill=lobato_set,x=log2FoldChange,geom='histogram')+facet_grid(lobato_set~.,scale='free_y')


# -

ggpubr::ggarrange(nrow=1,plotlist=list(
    ggMarginal(deneddata%>%filter(inner_q(Intensity,.99))%>%qplot(data=.,log10(1+Intensity),log10(1+`Intensity MG_1`),geom='blank')+  geom_point(),type='histogram'),
    deneddata%>%filter(inner_q(Intensity,.99))%>%qplot(data=.,1+`Intensity MG_2`,1+`Intensity MG_1`,log='xy')
    ))


library(orthologsBioMART)
orthos <- findOrthologsHsMm(from_filters = "hgnc_symbol",
  from_values = c("TP53","TERT"), 
  to_attributes = "external_gene_name")
print(orthos)

# ## Now compare the nedylation data to our proteomic fold changes



ned
