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

library(tidyverse)
library(data.table)
library(ggExtra)

library(here)
dened_protfile <- here('ext_data/dened_proteinGroups.txt')
if(!exists(dened_protfile)){
    dir.create(showWarn=F,here('ext_data'))
    download.file('ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/02/PXD022812/proteinGroups.txt',destfile = dened_protfile)
}
normalizePath(dened_protfile,mustWork=TRUE)

deneddata <- fread(dened_protfile)
selcols <- deneddata %>% colnames %>% str_subset(neg=TRUE,'counts|Razor|Unique')
selcols %>%paste0(col='\n')%>%message

# So it looks like we have 4 samples here, the MG_1 and 2, and NEDP1 and 2. I'll assume these are the proteotoxic and NEDD8 KO datasets respectively. Since there's no control group to compare to (it seems), we can just averag eintensities in our 

# +
options(repr.plot.width=21, repr.plot.height=7)
inner_q <- function(v,q=0.99){
    trim = (1 - q)/2
    (v >= quantile(v,trim)) & (v <= quantile(v,1-trim))
}
deneddata<-deneddata%>%mutate_at(vars(matches('Intensity')),as.numeric)
library(naniar)
library(ggExtra)

rmeanintens = deneddata%>%select(matches('Intensity.*_'))%>%as.matrix%>%rowMeans
qplot(deneddata$Intensity,rmeanintens)+ggtitle('the intensity column is just the mean of the other 4')


# -

# ## Harmonize gene identifiers
#
# We'll use a file downloaded from the biomart website to map our human gene names -> human ids -> mouse ids -> mouse names
#

deneddata$`Gene names`%>%head
fread(here('ext_data/human_mouse_orthologs_biomart.txt.gz'))%>%head(3)
deneddata%>%head

#

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
