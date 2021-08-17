
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
new_uniquedata<-filter(uniquedata[,76:107])
Genenames<- uniquedata$Gene.names
new_uniquedata$Gene.names<-Genenames
new_uniquedata<-new_uniquedata%>%
  relocate(Gene.names)
# new_uniquedata%>%across(vars(matches('_\\d')),median)
#also select the row with the highest median for each gene
gnmbestrows = new_uniquedata%>%
  ungroup%>%
  mutate(row=1:n())%>%
  pivot_longer(-one_of('row','Gene.names'))%>%
  group_by(Gene.names,row)%>%
  summarise_at(vars(value),median,na.rm=T)%>%
  group_by(Gene.names)%>%
  dplyr::slice(which.max(value))%>%
  .$row
new_uniquedata <- new_uniquedata[gnmbestrows,]
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

# View(new_uniquedata)
# View(uniquedata)
#Plot again as distribution of intensity (uniquepeptide_log2transform.png)

# new_uniquedata[-1]%>%
#   gather(variable, value)%>%
#   ggplot(aes(value))+
#   geom_histogram(bins=70) + 
#   facet_wrap(~variable, scales = 'free_y')

# #Plot as boxplot of distribution
# new_uniquedata[-1]%>%
#   gather(variable, value)%>%
#   ggplot(aes(x=value, y=variable))+
#   geom_boxplot(position="dodge2")

#remove batch effect using limma package where batch=plex1/2, batch2=replicates
#for replicates, a=3, b=2, c=1, d=5, e=4, f=ref)
#limma package requires BiocManager
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


gid2gname <- fread(here('pipeline/gid_2_gname.txt'),header=F)%>%
  set_colnames(c('gnm','g_id'))
# gid2gname <- setNames(gid2gname$gnm,gid2gname$g_id)

dds <- readRDS('data/dds.rds')
design(dds) <- 'group'
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
  c('MT','MT_Ctrl72'),
  c('MT','MT_A72')
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
  norm_med = norm_med,
)


normfunc=normfuncs[[1]]



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
      # compdata = tpm_grp_df%>%
      #   filter(sample_group==sampgroup)%>%
      #   inner_join(enframe(msdata,'gnm','ms'))%>%
      #   mutate_at(vars(tpm),log2)%>%
      #   # filter_at(vars(tpm,ms),is.finite)%>%
      #   identity
      compdata = rnaseqlfc%>%
        inner_join(enframe(mslfc,'gnm','ms_l2fc'),by='gnm')
      # tidy(cor.test(compdata$tpm,compdata$ms,use='complete'))
      cat('.')
      tidy(cor.test(compdata$log2FoldChange,compdata$ms_l2fc,use='complete'))
    #  qplot(data=compdata,log2FoldChange,ms_l2fc,alpha=I(0.5))  
    })
  })
})
})
stop()

foldchange_comps = map_df(.id='batchcor',batchcordsets,function(msdataset){
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
      # compdata = tpm_grp_df%>%
      #   filter(sample_group==sampgroup)%>%
      #   inner_join(enframe(msdata,'gnm','ms'))%>%
      #   mutate_at(vars(tpm),log2)%>%
      #   # filter_at(vars(tpm,ms),is.finite)%>%
      #   identity
      compdata = rnaseqlfc%>%
        inner_join(enframe(mslfc,'gnm','ms_l2fc'),by='gnm')
      # tidy(cor.test(compdata$tpm,compdata$ms,use='complete'))
      cat('.')
      tidy(cor.test(compdata$log2FoldChange,compdata$ms_l2fc,use='complete'))
    })
  })
})
})


foldchange_comps%>%select(batchcor,sample_pair,estimate,normfunc)%>%
  pivot_wider(values_from='estimate',names_from='sample_pair')



foldchange_comps

stop()


onebatchcor %>%
  gather(variable, value)%>%
  ggplot(aes(x=value, y=variable))+
  geom_boxplot(position="dodge2")+
  ggtitle("Single Batch corrected")

twobatchcor %>%
  gather(variable, value)%>%
  ggplot(aes(x=value, y=variable))+
  geom_boxplot(position="dodge2")+
  ggtitle("Double Batch corrected")

twobatchcor %>%
  gather(variable, value)%>%
  ggplot(aes(value))+
  geom_histogram(bins = 70)+
  ggtitle("Double Batch corrected")+
  facet_wrap(~variable, scales='free_x')

onebatchcor %>%
  gather(variable, value)%>%
  ggplot(aes(value))+
  geom_histogram(bins = 70)+
  ggtitle("Single Batch corrected")+
  facet_wrap(~variable, scales='free_x')




BNsingle <- bestNormalize(onebatchcor$Myoblast_3)


onebatchcor_copy <- onebatchcor
twobatchcor_copy <- twobatchcor

output <- apply(onebatchcor_copy,2, bestNormalize)
sink("bestNormalizeOutput.txt")
print(output)
sink()


#Repeat using twobatch correction data
output2 <- apply(twobatchcor_copy,2, bestNormalize)
sink("bestNormalizeOutput2.txt")
print(output)
sink()

#Perform OrderNorm on all columns which seems to fit most samples as per bestNormalizeOutput2.txt
test <- orderNorm(twobatchcor$Myoblast_3)
test
libary(MASS)
MASS::truehist(test$x.t, main="OrderNorm transformation", nbins=30)
myob3_ON<-test$x.t
myob3_ON<- as.data.frame(myob3_ON)
ggplot(as.data.frame(test$x.t), aes(test$x.t))+
  geom_boxplot()+
  ggtitle("OrderNorm transformation")

myob3_ON<-fortify(myob3_ON)

#create a list of doubles containing the $x.t (transformed data) in output2
normdata <- lapply(output2, '[[', 1)
# View(normdata)

df1 <- data.frame(matrix("", ncol = 0, nrow = length(twobatchcor$Myoblast_3)))

for(i in seq(1,length(colnames(twobatchcor)))) {
  str01 <- paste(colnames(twobatchcor[i]), sep = "")
  templst <- orderNorm(twobatchcor[,i])
  df1[,i] <- as.data.frame(templst$x.t)
  colnames(df1)[i] <- str01
}



df1$Gene.names <- Gene.names
df1<-df1%>%
  relocate(Gene.names)
head(df1)

#Plot distribution of orderNorm values to check normalization
df1[-1]%>%
  gather(variable, value)%>%
  ggplot(aes(value))+
  geom_histogram(bins=70) + 
  facet_wrap(~variable, scales = 'free_y')

df1_copy <-select(df1,-(Gene.names))

df1_copy %>%
  gather(variable, value)%>%
  ggplot(mapping=aes(x=value))+
  geom_histogram(mapping=aes(colour=variable), bins=100)

df1_long <- df1 %>%
  pivot_longer(!Gene.names, names_to = "Treatment", values_to = "Intensity")

head(df1_long) 

df1_long%>%
  filter(Gene.names=="Dcaf6")%>%
  ggplot(aes(x=Treatment, y=Intensity))+
  geom_point()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

#Selecting data by gene name and by treatment
# df1_long%>%
#   filter(Gene.names=="Trim63")%>%
#   filter(str_detect(Treatment, "Myoblast"))



trim63_data <- df1_long%>%
  filter(Gene.names=="Nae1")%>%
  filter(str_detect(Treatment, "^Myoblast"))

test <- df1 %>%
  gather(replicate, value, -1) %>%
  #filter(!duplicated(replicate)) %>%
  filter(Gene.names == "Nae1") %>%
  filter((replicate != "Ref1") & (replicate != "Ref2")) %>%
  arrange(replicate) %>%
  mutate(replicate = substr(replicate,1,nchar(replicate)-2))

# head(test)
# test %>%
#     summarise(mean=mean(value))

sample <- c("Myoblast", "Myotube", "Atrophy_Ctrl_24h", "Atrophy_Ctrl_72h", "Atrophy_24h", "Atrophy_72h")
plot<-test %>%
  group_by(replicate)%>%
  ggplot(aes(x=factor(replicate, level=sample), y=value),colour=replicate)+
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=8, size=8, color="red", fill="red")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  ggtitle("Cul1")+
  labs(x="Treatment",y="Normalized Intensity")

png(filename="Cul1_box.png", width=800, height = 800, units="px", pointsize=12)
plot
dev.off()


GOIs <- c("Trim63", "Myh1")

# df1 %>%
#   gather(replicate, value, -Gene.names) %>%
#   mutate(replicate = 
#            ifelse((replicate != "Ref1") & (replicate != "Ref2"),
#                   substr(replicate,1,nchar(replicate)-2), replicate)) %>%
#   filter(Gene.names == GOIs[1]) %>%
#   ggplot(aes(x = replicate, y = value)) + 
#   geom_boxplot(middle=mean()) +
#   labs(title=GOIs[1])+
#   theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

