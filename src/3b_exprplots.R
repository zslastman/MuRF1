library(tidyverse)
library(magrittr)
library(here)
library(DESeq2)
library(data.table)

samporder =c("MB","MT","MT_Ctrl24","MT_Ctrl72","MT_A24","MT_A72")
gnm_gid = fread('data/gid_2_gname.txt',header=F)
tx_countdata <- readRDS(here('data/tx_countdata.rds'))
# genesofinterest = read_tsv('ext_data/2000929_TranscriptomeCandList.txt')

genesofinterest = fread('ext_data/genes_of_interest_manual.tsv',header=F)%>%
  set_colnames(c('gene_name','p_ids'))
newgns = 'ext_data/Expr_gene_Jan2021.txt'%>%fread(header=F)%>%set_colnames('gene_name')%>%
  mutate(p_ids = '_')%>%
  mutate(gene_name = str_to_sentence(gene_name))

stopifnot(all(newgns$gene_name%>%is_in(gnm_gid$V1)))
# genesofinterest %<>%bind_rows(newgns)
genesofinterest <- (newgns)
normcounts<-readRDS(here('data/tx_countdata.rds'))
dds <- readRDS(here('data/dds.rds'))
dds <- DESeq(dds)

resultslist <- readRDS(here('data/resultslist.rds'))

modelmat = model.matrix(data=colData(dds),design(dds))

i = 1
pred_conf_df = map_df(1:nrow(unique(modelmat)),function(i){
  #    
  contrnm = modelmat%>%unique%>%rownames%>%.[i]
  contrnm = contrnm%>%str_replace('_1','') 
  #
  contrv  = modelmat%>%unique%>%.[i,]
  #
  res = results(dds,contrast = contrv)
  #
  lfcdf = tibble(g_id = rownames(res),mean=res$log2FoldChange,lfc_se = res$lfcSE)%>%
    mutate(sample = contrnm)%>%
    left_join(gnm_gid%>%set_colnames(c('gene_name','g_id')),by='g_id')
  lfcdf
})


dir.create( here(paste0('plots/geneexprplots/')))
pdf = grDevices::pdf
for(igene in genesofinterest$gene_name){
# for(igene in 'Cul4a'){
  pids = genesofinterest$p_ids[match(igene,genesofinterest$gene_name)]%T>%{stopifnot(!is.na(.))}
  gid =  gnm_gid$V2[match(igene,gnm_gid$V1)] %T>%{stopifnot(!is.na(.))}
  #
  expr = tx_countdata$abundance[gid,]%>%enframe('sample','TPM')%>%mutate(gene_name=igene)
  expr$treatment = expr$sample%>%str_extract('Ctrl|A(?=\\d)')%>%replace_na('N')
  expr$sample = expr$sample%>%str_replace('_\\d+$','')
  expr%<>%arrange(match(sample,samporder))
  expr$sample%<>%as_factor
  expr$treatment%<>%as_factor
  # 
  errorbar_dat = expr%>%group_by(gene_name,treatment,sample)%>%summarise(TPM=mean(TPM))%>%
    left_join(pred_conf_df,by=c('sample','gene_name'))%>%
    mutate(
      lower= 2^(log2(TPM) -1.96*lfc_se),
      upper= 2^(log2(TPM) + 1.96*lfc_se)
    )
  #now plot
  plotfile<- here(paste0('plots/geneexprplots/',igene,'.',pids,'exprplot','.pdf'))
  pdf(plotfile,w=7,h=4)
  p=expr%>%
    ggplot(.,aes(x=sample,y=TPM,color=treatment))+
    geom_point(position=position_jitter(width=0.2))+
    geom_linerange(data=errorbar_dat,aes(ymin=lower,ymax=upper))+
    # geom_(data=errorbar_dat,aes(ymin=lower,ymax=upper))+
    # scale_color_discrete(name='colorname',colorvals)+
    scale_x_discrete(paste0('Sample'))+
    scale_y_continuous(paste0('TPM'))+
    ggtitle(paste0(igene,'\n',pids))+
    theme_bw()
  print(p)
  dev.off()
  message(normalizePath(plotfile))
}



