
library(tidyverse)
library(here)
library(magrittr)
library(ggplot2)
library(BiocManager)
library(bestNormalize)
library(dplyr)
library(limma)
library(DESeq2)
# BiocManager::install("version=3.12")

rename <- dplyr::rename
slice <- dplyr::slice
# setwd("C:/Users/snakaga/MDC/Labbook/MassSpec/Differential time course and atrophy/normalization_struggles")
data<-read.delim("ext_data/proteinGroups_withNewPlex1.txt", sep="\t", header=TRUE)

colnames(data)
#reverse matches and contaminants are entered as + so remove anything not an empty cell
#check dplyr attached, remove + entries in Reverse and Potential.contaminant columns
cleandata <- data %>%
  filter(Reverse == '')%>%
  filter(Potential.contaminant == '')
# source("", chdir = TRUE)
#rename columns to match channel and plex assignment
cleandata <- rename(cleandata, Myoblast_3= Reporter.intensity.corrected.1.Plex1,
    Myotube_2=Reporter.intensity.corrected.2.Plex1,
  Atrophy_Ctrl_72h_2=Reporter.intensity.corrected.3.Plex1,
  Atrophy_Ctrl_24h_2=Reporter.intensity.corrected.4.Plex1,
  Atrophy_72h_3=Reporter.intensity.corrected.5.Plex1,
  Myoblast_1=Reporter.intensity.corrected.6.Plex1,
  Atrophy_24h_5=Reporter.intensity.corrected.7.Plex1,
  Atrophy_72h_5=Reporter.intensity.corrected.8.Plex1,
  Atrophy_72h_1=Reporter.intensity.corrected.9.Plex1,
  Myotube_4=Reporter.intensity.corrected.10.Plex1,
  Myoblast_5=Reporter.intensity.corrected.11.Plex1,
  Atrophy_24h_3=Reporter.intensity.corrected.12.Plex1,
  Atrophy_Ctrl_24h_4=Reporter.intensity.corrected.13.Plex1,
  Atrophy_Ctrl_72h_4=Reporter.intensity.corrected.14.Plex1,
  Atrophy_24h_1=Reporter.intensity.corrected.15.Plex1,
  ef1=Reporter.intensity.corrected.16.Plex1)
cleandata <- rename(cleandata,
    Myotube_5=Reporter.intensity.corrected.1.Plex2,
  Atrophy_Ctrl_72h_5=Reporter.intensity.corrected.2.Plex2,
  Myoblast_2=Reporter.intensity.corrected.3.Plex2,
  Atrophy_72h_4=Reporter.intensity.corrected.4.Plex2,
  Myotube_3=Reporter.intensity.corrected.5.Plex2,
  Atrophy_Ctrl_24h_3=Reporter.intensity.corrected.6.Plex2,
  Atrophy_Ctrl_24h_5=Reporter.intensity.corrected.7.Plex2,
  Myoblast_4=Reporter.intensity.corrected.8.Plex2,
  Atrophy_Ctrl_72h_3=Reporter.intensity.corrected.9.Plex2,
  Atrophy_24h_2=Reporter.intensity.corrected.10.Plex2,
  Myotube_1=Reporter.intensity.corrected.11.Plex2,
  Atrophy_Ctrl_24h_1=Reporter.intensity.corrected.12.Plex2,
  Atrophy_24h_4=Reporter.intensity.corrected.13.Plex2,
  Atrophy_Ctrl_72h_1=Reporter.intensity.corrected.14.Plex2,
  Atrophy_72h_2=Reporter.intensity.corrected.15.Plex2,
  Ref2=Reporter.intensity.corrected.16.Plex2)
colnames(cleandata)

#check all sample groups have 5 replicates each
length(grep(x=colnames(cleandata),pattern="Myoblast"))
length(grep(x=colnames(cleandata), pattern="Myotube"))
length(grep(x=colnames(cleandata), pattern="Atrophy_Ctrl_24h"))
length(grep(x=colnames(cleandata), pattern="Atrophy_Ctrl_72h"))
length(grep(x=colnames(cleandata), pattern="Atrophy_24h"))
length(grep(x=colnames(cleandata), pattern="Atrophy_72h"))

# View(cleandata)

newdata <- filter(cleandata[,76:107])
newdata$Gene.names <- cleandata$Gene.names
ncol(newdata)
colnames(newdata)

newdata%>%
  relocate(Gene.names)

# View(cleandata)

#Create uniquedata filtering by unique.peptides >=2 in both
uniquedata<- cleandata %>%
  filter(Unique.peptides.Plex1>=2|Unique.peptides.Plex2>=2)
#Create notuniquedata filtering by peptides including non-unique peptides
notuniquedata <- cleandata %>%
  filter(Peptides.Plex1>=2|Peptides.Plex2>=2)

#Select only treatment group columns and add genenames column
udata_sel <-filter(uniquedata[,76:107])
Genenames<- uniquedata$Gene.names
udata_sel $Gene.names<-Genenames
udata_sel <-udata_sel %>%
  relocate(Gene.names)
# udata_sel %>%across(vars(matches('_\\d')),median)
#also select the row with the highest median for each gene
gnmbestrows = udata_sel %>%
  ungroup%>%
  mutate(row=1:n())%>%
  pivot_longer(-one_of('row','Gene.names'))%>%
  group_by(Gene.names,row)%>%
  summarise_at(vars(value),median,na.rm=T)%>%
  group_by(Gene.names)%>%
  dplyr::slice(which.max(value))%>%
  .$row
new_uniquedata  <- udata_sel [gnmbestrows,]
new_uniquedata$Gene.names%<>%as.character
#
stopifnot(!new_uniquedata$Gene.names%>%duplicated%>%any)
# #Plot distribution of intensity using script line 34~41
# new_uniquedata[-33]%>%
#   gather(variable, value)%>%
#   ggplot(aes(value))+
#   geom_histogram(bins = 70) + 
#   xlim(-1,5e+06)+
#   facet_wrap(~variable, scales = 'free_x')
#log2 transformation to prepare for limma package
colnames(new_uniquedata)
new_uniquedata[,-1] <-log(new_uniquedata[,-1],2) 
#remove -Inf entries i.e. those w/ 0.00 prior to log2 transform
new_uniquedata<-new_uniquedata%>%
  filter_at(vars(2:33), all_vars(!is.infinite(.)))


colnames(new_uniquedata[,-1])==colnames(udata_sel[,-1])


batch <- c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B")
batch2 <- c("A", "B", "B", "B", "A", "C", "D","D","C", "E", "D", "A", "E", "E", "C", "F", "E","E","B", "E", "A","A", "D", "E", "A", "B", "C","C","E","C", "B", "F") 
is_myob <- colnames(new_uniquedata[,-1])%>%str_detect('Myoblast')
# batch = paste0(batch,ifelse(is_myob,'MB','MT'))
# batch2 = paste0(batch2,ifelse(is_myob,'MB','MT'))
new_uniquedata%>%colnames
#check for presence of na in dataframe
#checking for both na and -Inf/Inf would be apply(new_uniquedata, 2, function(x) any(is.na(x)|is.infinite(x)))
indx <-apply(new_uniquedata, 2, function(x) any(is.na(x)))
indx
nobatchcor <- new_uniquedata[,-1]
twobatchcor<- new_uniquedata[-1]%>%
  removeBatchEffect(batch, batch2)
onebatchcor<- new_uniquedata[-1]%>%
  removeBatchEffect(batch)
nobatchcor<- as.data.frame(nobatchcor)%>%
  set_rownames(new_uniquedata[[1]])
onebatchcor<- as.data.frame(onebatchcor)%>%
  set_rownames(new_uniquedata[[1]])
twobatchcor <- as.data.frame(twobatchcor)%>%
  set_rownames(new_uniquedata[,1])


library(data.table)
library(here)
library(magrittr)
library(broom)


gid2gname <- read_tsv(here('pipeline/gid_2_gname.txt'),col_names=F)%>%
  set_colnames(c('gnm','g_id'))
# gid2gname <- setNames(gid2gname$gnm,gid2gname$g_id)

dds <- readRDS('data/dds.rds')
design(dds) <- as.formula('~group')
dds <- DESeq(dds)
resnames <- c("Intercept", "group_MT_vs_MB", "group_MT_A24_vs_MB", "group_MT_A72_vs_MB", 
"group_MT_Ctrl24_vs_MB", "group_MT_Ctrl72_vs_MB")
stopifnot(resultsNames(dds)==resnames)
grpres = lapply(resnames,function(rs) results(dds,contrast=list(rs)))
grpnames = resnames%>%str_replace('Intercept','MB')%>%str_replace('_vs_MB','')%>%str_replace('group_','')
grpres%>%setNames(grpnames)%>%map_df(~as.data.frame(.)%>%rownames_to_column('g_id')%>%select(g_id,))

tx_countdata <- 'data/tx_countdata.rds'%>%readRDS 
tpmdf = tx_countdata%>%.$abundance%>%as.data.frame%>%rownames_to_column('g_id')%>%
  pivot_longer(-g_id,names_to='sample',values_to='tpm')%>%
  left_join(gid2gname,by='g_id')
tpmdf$sample_group = tpmdf$sample%>%str_replace('_\\d+$','')
tpm_grp_df<- tpmdf%>%group_by(sample_group,gnm)%>%summarise_at(vars(tpm),mean)

sampgroups <- tpmdf$sampgroup%>%unique
sampgroup <- sampgroups[[1]]
c("MB", "MT", "MT_A24", "MT_A72", "MT_Ctrl24", "MT_Ctrl72")

{
mscolnames_replace <- function(x){
  x%>%
     str_replace('Myoblast','MB')%>%
    str_replace('Myotube','MT')%>%
    str_replace('Atrophy_Ctrl_(\\d+)h','MT_Ctrl\\1')%>%
    str_replace('Atrophy_(\\d+)h','MT_A\\1')
}
mscolnames_replace
}


sampgroup_pairs = list(
  c('MB','MT'),
  c('MT','MT_Ctrl24'),
  c('MT','MT_A24'),
  c('MT_Ctrl24','MT_Ctrl72'),
  c('MT_Ctrl24','MT_A24'),
  c('MT_Ctrl24','MT_Ctrl72'),
  c('MT_Ctrl72','MT_A72'),
  c('MT_A24','MT_A72')
)%>%setNames(.,map_chr(.,paste0,collapse='_'))

norm_MAD <- function(x) {
  (x - median(x, na.rm=T))/mad(x, na.rm=T)
}

norm_med <- function(x) {
  x - median(x, na.rm=T)
}

batchcordsets = list(nobatchcor,onebatchcor,twobatchcor)%>%
  setNames(c('nobatchcor','onebatchcor','twobatchcor'))
# comparisontypes = c('intrasample','lfc2')
comparisontypes = c('lfc2')%>%setNames(.,.)
msdataset=batchcordsets[[2]]
comparisontype=comparisontypes[[1]]
samppair = sampgroup_pairs[[1]]
normfuncs = c(
  'none'=identity,
  'bestNormalize'=function(x)bestNormalize::orderNorm(x)$x.t,
  norm_MAD = norm_MAD,
  norm_med = norm_med
)

samppair=sampgroup_pairs[[1]]
contrnm=names(sampgroup_pairs)[[1]]
normfunc=identity
library(DESeq2)
foldchange_comps = map_df(.id='batchcor',batchcordsets['onebatchcor'],function(msdataset){
  map_df(.id='comptype',comparisontypes,function(comparisontype){
  map_df(.id='normfunc',normfuncs,function(normfunc){
    # lapply(sampgroup,function(stage){
    map_df(.id='sample_pair',sampgroup_pairs,function(samppair){
      # msdataset%<>%head
      samp1=samppair[[1]]
      samp2=samppair[[2]]
      msdatacols1 = msdataset%>%colnames%>%
        mscolnames_replace%>%
        str_detect(paste0(samp1,'_\\d'))
      msdatacols2 = msdataset%>%colnames%>%
        mscolnames_replace%>%
        str_detect(paste0(samp2,'_\\d'))
      stopifnot(any(msdatacols1))
      stopifnot(any(msdatacols2))
      msdata1 = msdataset[,msdatacols1]%>%rowMeans(na.rm=T)%>%
        setNames(rownames(msdataset))
      msdata2 = msdataset[,msdatacols2]%>%rowMeans(na.rm=T)%>%
        setNames(rownames(msdataset))
      mslfc = (normfunc(msdata2) - normfunc(msdata1))
      #
      rnaseqlfc = results(dds,c('group',samp2,samp1))
      rnaseqlfc= rnaseqlfc%>%as.data.frame%>%rownames_to_column('g_id')%>%
        left_join(gid2gname,by='g_id')%>%
        filter(padj<0.05)%>%
        select(gnm,log2FoldChange)
      compdata = rnaseqlfc%>%
        inner_join(enframe(mslfc,'gnm','ms_l2fc'),by='gnm')
      # tidy(cor.test(compdata$tpm,compdata$ms,use='complete'))
      cat('.')
      tidy(cor.test(compdata$log2FoldChange,compdata$ms_l2fc,use='complete'))
     # qplot(data=compdata,log2FoldChange,ms_l2fc,alpha=I(0.5))  
    })
  })
})
})
dir.create('tables')

foldchange_comps%>%write_tsv('tables/foldchange_correlations.tsv')

imap(normfuncs[c('none','norm_MAD','norm_med')],function(normfunc,normfuncname){
udata_sel%>%
  mutate_at(vars(-Gene.names),log2)%>%
  {set_rownames(as.matrix(.[,-1]),.[,1])}%>%
  removeBatchEffect(batch)%>%
  normfunc%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  write_tsv(str_interp('tables/tmt_ms_${normfuncname}.tsv'))
})