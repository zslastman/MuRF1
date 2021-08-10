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
# - print fold change plots
# - output some tables
#

library(tidyverse);
library(data.table);
library(magrittr);

# +
library(here)
#readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3)
mg132_gnms = suppressMessages(readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3))[[2]]%>%.[!is.na(.)]
nedd8ko_gnms = suppressMessages(readxl::read_xlsx('../ext_data/lobatogil_etal_2021_s1.xlsx',sheet=3))[[5]]

message(str_interp('${length(mg132_gnms)} mg132 gene names'))
message(str_interp('${length(nedd8ko_gnms)} nedd8 ko gene names'))



lobato_gnms = bind_rows(.id='lobato_set', mg132=tibble(gnm=mg132_gnms), nedd8ko = tibble(gnm=nedd8ko_gnms))
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

# # Load mRNA fold changes and intersect with lobato IDs

resultslist <- readRDS(here('data/resultslist.rds'))
mrna_gids = resultslist[[1]]$gene_id
message('ow many of these genes have at least one ortholog appearing our annotation?')
lobato_gnm_orth%>%group_by(gnm)%>%summarise(inmrna_dat = any(mgid %in% mrna_gids))%>%.$inmrna_dat%>%table


# ## Now print plots comparing our mRNA fold changes for the two sets to the general distribution

# +
options(jupyter.plot_mimetypes = "image/svg+xml") 

options(repr.plot.width=21, repr.plot.height=7)
inner_q <- function(v,q=0.99){
    trim = (1 - q)/2
    (v >= quantile(v,trim)) & (v <= quantile(v,1-trim))
}

squish <- function(v,q=0.99){
    trim = (1 - q)/2
    lower = quantile(v,trim)
    upper = quantile(v,1-trim)
    v%>%pmax(lower)%>%pmin(upper)
}


contrasts = names(resultslist)
contrast = contrasts[[1]]

get_lobatoset_histogram <- function(contrast){
    lobato_lfc_df = resultslist[[contrast]]%>%
        left_join(lobato_gnm_orth,by=c('gene_id'='mgid'))%>%
        distinct(lobato_set,gnm,log2FoldChange,padj)
    
    # lobato_lfc_df%>%sample_n(4)

    # lobato_lfc_df%>%
    #     filter(!is.na(log2FoldChange))%>%
    #     filter(!is.na(lobato_set)|inner_q(log2FoldChange,0.99))%>%
    #     qplot(data=.,fill=lobato_set,x=log2FoldChange,geom='histogram')+facet_grid(lobato_set~.,scale='free_y')+
    #     ggtitle(str_interp('distribution of Neddylation classes for ${contrast}'))

    
    ggdf=lobato_lfc_df%>%
        filter(padj<0.05)%>%
        mutate(log2FoldChange = squish(log2FoldChange))
    ggdf$lobato_set%<>%replace_na('other')
    lfclist = split(ggdf$log2FoldChange,ggdf$lobato_set)
    library(broom)
    mg132less = wilcox.test(alternative = 'less',lfclist[['mg132']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    mg132more = wilcox.test(alternative = 'greater',lfclist[['mg132']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    nedd8less = wilcox.test(alternative = 'less',lfclist[['nedd8ko']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    nedd8more = wilcox.test(alternative = 'greater',lfclist[['nedd8ko']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}

    

    ggdf%>%
        qplot(data=.,fill=lobato_set,x=log2FoldChange,geom='blank')+facet_grid(lobato_set~.,scale='free_y')+
        geom_histogram(binwidth=0.5)+
        scale_x_continuous(breaks=function(l)seq(floor(l[1]),ceiling(l[2])))+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=mg132less,x=-Inf,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=mg132more,x=Inf,y=Inf),hjust=1,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=nedd8less,x=-Inf,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=nedd8more,x=Inf,y=Inf),hjust=1,vjust=1)+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=length(lfclist$mg132),x=0,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=length(lfclist$nedd8ko),x=0,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='other'),aes(label=length(lfclist[['other']]),x=0,y=Inf),hjust=0,vjust=1)+
        ggtitle(str_interp('distribution of mRNA log fold change for nedlation classes, \n contrast: ${contrast}'))
    #     coord_cartesian(xlim=c(-10,10))
}
plots = lapply(contrasts,get_lobatoset_histogram)
ggpubr::ggarrange(plotlist=plots,ncol=4)
# -

# ## Now compare the nedylation data to our proteomic fold changes

# +
tmtdata = fread(here('ext_data/Results_Proteome_20210115.txt'))

gid2gname <- fread(here('pipeline/gid_2_gname.txt'),header=F)%>%set_colnames(c('gnm','g_id'))
mtmtdata = tmtdata%>%dplyr::rename('gnm':=id)%>%filter(gnm %in% gid2gname$gnm)



tmtcols <- c(
  "logFC.diff_atrophy_8.over.diff_atrophy_10",
  "logFC.diff_Ctr_10.over.diff_atrophy_10",
  "logFC.diff_Ctr_7.over.diff_atrophy_10",
  "logFC.diff_Ctr_8.over.diff_atrophy_10",
  "logFC.undiff_Ctr_0.over.diff_atrophy_10",
  "logFC.diff_Ctr_10.over.diff_atrophy_8",
  "logFC.diff_Ctr_7.over.diff_atrophy_8",
  "logFC.diff_Ctr_8.over.diff_atrophy_8",
  "logFC.undiff_Ctr_0.over.diff_atrophy_8",
 "logFC.diff_Ctr_7.over.diff_Ctr_10",
 "logFC.diff_Ctr_8.over.diff_Ctr_10",
 "logFC.undiff_Ctr_0.over.diff_Ctr_10",
 "logFC.diff_Ctr_8.over.diff_Ctr_7",
 "logFC.undiff_Ctr_0.over.diff_Ctr_7",
 "logFC.undiff_Ctr_0.over.diff_Ctr_8"
)%>%setNames(.,.)


# +
tmtcol = tmtcols[[1]]
tmtpcol=tmtcol%>%str_replace('logFC','adj.P.Val')

get_lset_prot_histogram <- function(tmtcol){
    tmtpcol=tmtcol%>%str_replace('logFC','adj.P.Val')
    plobato_lfc_df = mtmtdata%>%select(gnm,log2FoldChange=!!tmtcol,padj = !!tmtpcol)%>%
        left_join(lobato_gnm_orth%>%dplyr::rename('mgnm':=`Mouse gene name`),by=c('gnm'='mgnm'))%>%
        distinct(lobato_set,gnm,log2FoldChange,padj)
    

    # lobato_lfc_df%>%sample_n(4)

    # lobato_lfc_df%>%
    #     filter(!is.na(log2FoldChange))%>%
    #     filter(!is.na(lobato_set)|inner_q(log2FoldChange,0.99))%>%
    #     qplot(data=.,fill=lobato_set,x=log2FoldChange,geom='histogram')+facet_grid(lobato_set~.,scale='free_y')+
    #     ggtitle(str_interp('distribution of Neddylation classes for ${contrast}'))

    
    ggdf=plobato_lfc_df%>%
        filter(padj<0.05)%>%
        mutate(log2FoldChange = squish(log2FoldChange))
    ggdf$lobato_set%<>%replace_na('other')
    lfclist = split(ggdf$log2FoldChange,ggdf$lobato_set)
    library(broom)
    mg132less = wilcox.test(alternative = 'less',lfclist[['mg132']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    mg132more = wilcox.test(alternative = 'greater',lfclist[['mg132']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    nedd8less = wilcox.test(alternative = 'less',lfclist[['nedd8ko']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}
    nedd8more = wilcox.test(alternative = 'greater',lfclist[['nedd8ko']],lfclist[['other']])%>%tidy%>%{formatC(.$p.value, format = "e", digits = 2)}

    

    plot = ggdf%>%
        qplot(data=.,fill=lobato_set,x=log2FoldChange,geom='blank')+facet_grid(lobato_set~.,scale='free_y')+
        geom_histogram(binwidth=0.5)+
        scale_x_continuous(breaks=function(l)seq(floor(l[1]),ceiling(l[2])))+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=mg132less,x=-Inf,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=mg132more,x=Inf,y=Inf),hjust=1,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=nedd8less,x=-Inf,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=nedd8more,x=Inf,y=Inf),hjust=1,vjust=1)+
        geom_text(data=tibble(lobato_set='mg132'),aes(label=length(lfclist$mg132),x=0,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='nedd8ko'),aes(label=length(lfclist$nedd8ko),x=0,y=Inf),hjust=0,vjust=1)+
        geom_text(data=tibble(lobato_set='other'),aes(label=length(lfclist[['other']]),x=0,y=Inf),hjust=0,vjust=1)+
        ggtitle(str_interp('distribution of mRNA log fold change for nedlation classes, \n contrast: ${tmtcol}'))
    #     coord_cartesian(xlim=c(-10,10))
    list(plot,ggdf)
}
plots = lapply(tmtcols,get_lset_prot_histogram)
ggdfs = plots%>%map(2)
plots = plots%>%map(1)
ggpubr::ggarrange(plotlist=plots,ncol=4)
                           
# -


# ## Finally, export table con

# ggdfs%>%bind_rows(.id='MS_contrast')%>%filter(lobato_set!='other')%>%spread(MS_contrast)
ggdfs%>%bind_rows(.id='MS_contrast')%>%filter(lobato_set!='other')%>%write_tsv('../tables/lobato_tmt.tsv')
# + active=""
# ## we want to know what percentage of the lobato sets are up/down/either regulated at each timepoint
# +

# lobato_gnm_orth%>%filter('%>%distinct(gnm)%>%left_join()

iset = lobato_gnm_orth$lobato_set%>%unique%>%head(1)

isetgnms = lobato_gnm_orth%>%filter(lobato_set==iset)%>%select(gnm=`Mouse gene name`)%>%distinct%>%filter(gnm!='')%>%.[[1]]

isets = lobato_gnm_orth$lobato_set%>%unique%>%head(2)
updownpcdf = map_df(.id='lobato_set',isets%>%setNames(isets),function(iset){
 ggdfs%>%bind_rows(.id='MS_contrast')%>%filter(gnm %in% isetgnms)%>%
    mutate(dir=ifelse(log2FoldChange>0,'up','down'))%>%
    group_by(MS_contrast,dir)%>%tally%>%
    mutate(pc_set = n / length(isetgnms))
})

updownpcdf%>%write_tsv('../tables/lobato_pc_reg.tsv')

# lobato_gnm_orth%>%select(lobato_set,gnm=`Mouse gene name`)%>%distinct%>%left_join(ggdfs%>%bind_rows(.id='MS_contrast'))%>%
#     mutate(dir=ifelse(log2FoldChange>0,'up','down'))%>%select(lobato_set,gnm,MS_contrast,dir)%>%mutate(tmp=TRUE)%>%spread(dir,tmp)%>%
#     select(-`<NA>`)%>%
#     group_by(MS_contrast)
# -

updownpcdf$MS_contrast%<>%str_replace('logFC.(un)?diff_','')%>%str_replace('.over.diff.','.vs.\n')
updownpcdf %>% ggplot(data=.,aes(y=pc_set,x=dir,fill=dir))+stat_identity(geom='bar')+facet_grid(lobato_set~MS_contrast)


