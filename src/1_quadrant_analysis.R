library(tidyverse)
library(magrittr)
library(limma)
library(here)
library(Biobase)
library(proDA)
library(DESeq2)

slice <- dplyr::slice

gid2gnm <- read_tsv('pipeline/gid_2_gname.txt',col_names=c('gene_name','gene_id'))
gid2gnm <- setNames(gid2gnm$gene_name,gid2gnm$gene_id)

{
#
sampgroup_pairs = list(
  c('MB','MT'),
  c('MT','MT_Ctrl24'),
  c('MT','MT_A24'),
  c('MT_Ctrl24','MT_Ctrl72'),
  c('MT_A24','MT_A72'),
  c('MT_Ctrl24','MT_A24'),
  c('MT_Ctrl72','MT_A72')
)%>%setNames(.,map_chr(.,paste0,collapse='_'))
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
#
ms_data <- read_tsv('tables/tmt_ms_norm_med.tsv')%>%
	select(-matches('Ref'))%>%
	select(-matches('ef1'))
colnames(ms_data)%<>%mscolnames_replace
#
colnames(ms_data)[-1]%<>%paste0('MS_',.)
stopifnot(colnames(ms_data)==
	c("gene_name", "MS_MB_3", "MS_MT_2", "MS_MT_Ctrl72_2", "MS_MT_Ctrl24_2", 
"MS_MT_A72_3", "MS_MB_1", "MS_MT_A24_5", "MS_MT_A72_5", "MS_MT_A72_1", "MS_MT_4", 
"MS_MB_5", "MS_MT_A24_3", "MS_MT_Ctrl24_4", "MS_MT_Ctrl72_4", "MS_MT_A24_1", 
 "MS_MT_5", "MS_MT_Ctrl72_5", "MS_MB_2", "MS_MT_A72_4", "MS_MT_3",
"MS_MT_Ctrl24_3", "MS_MT_Ctrl24_5", "MS_MB_4", "MS_MT_Ctrl72_3", "MS_MT_A24_2",
"MS_MT_1", "MS_MT_Ctrl24_1", "MS_MT_A24_4", "MS_MT_Ctrl72_1",
"MS_MT_A72_2")
)
library(txtplot)
ms_data[[2]]%>%.[is.finite(.)]%>%txtdensity
msmedians = ms_data[,-1]%>%apply(2,median)
#median norm
ms_data[,-1] <- ms_data[,-1]%>%sweep(2,F='-',STATS=msmedians)
}


################################################################################
########## Model RNAseq with Limma
################################################################################
#
{
	gid2gnm <- read_tsv('pipeline/gid_2_gname.txt',col_names=c('gene_name','gene_id'))
	gid2gnm <- setNames(gid2gnm$gene_name,gid2gnm$gene_id)
	tx_countdata <- 'data/tx_countdata.rds'%>%readRDS
	colnames(tx_countdata$counts)%<>%paste0('rnaseq_',.)
	rownames(tx_countdata$counts)<-gid2gnm[rownames(tx_countdata$counts)]
	#
	tx_countdata$counts%>%rowMedians%>%add(1)%>%log2%>%txtdensity
	lowcounts <- tx_countdata$counts%>%rowMedians%>%`<`(1)
	tx_countdata$counts<-tx_countdata$counts[!lowcounts,]
	countdesign <- tibble(
		sample_id = tx_countdata$counts%>%colnames,
		group = tx_countdata$counts%>%colnames%>%str_replace('_\\d$','')%>%
			str_replace('rnaseq_',''),
		assay='rnaseq'
	)
	countvoom <- voom(tx_countdata$counts, model.matrix(data=countdesign,~group))
	countfit <- lmFit(countvoom)%>%eBayes
	#
	countvoom <- voom(tx_countdata$counts,
		model.matrix(data=countdesign,~group))
	#
	#
	#
	msdesign <- tibble(
		sample_id = ms_data[,-1]%>%colnames,
		group = ms_data[,-1]%>%colnames%>%str_replace('_\\d$','')%>%
			str_replace('MS_',''),
		assay='MS'
	)
	grplevels <- c("MB", "MT", "MT_Ctrl24", "MT_Ctrl72","MT_A24", "MT_A72")
	msdesign$group%<>%factor(levels=grplevels)
	#
	msvoom <- ms_data%>%
		select(-matches('Ref'))%>%
		{set_rownames(as.matrix(.[,-1]),.[[1]])}%>%
		{exp(.)}%>%
		voom(design=model.matrix(data=msdesign,~group))
	msvoom$weights%<>%set_colnames(colnames(msvoom$E))
	msvoom$design%<>%set_rownames(colnames(msvoom$E))
	msvoom$targets%<>%set_rownames(colnames(msvoom$E))
	msfit <- lmFit(msvoom)%>%eBayes
	#
	sharedgnms <- intersect(rownames(countvoom$E),rownames(msvoom$E))
	countvoom%>%names
	rownames(msvoom$weights)<-rownames(msvoom$E)
	rownames(countvoom$weights)<-rownames(countvoom$E)
	alldesign <- bind_rows(countdesign,msdesign)
	grplevels <- c("MB", "MT", "MT_Ctrl24", "MT_Ctrl72","MT_A24", "MT_A72")
	alldesign$group%<>%factor(levels=grplevels)
	alldesign$assay%<>%factor(levels=c('rnaseq','MS'))
	#
	allvoom<-new('EList')
	allvoom$E<-cbind(countvoom$E[sharedgnms,],msvoom$E[sharedgnms,])
	allvoom$weights<-cbind(countvoom$weights[sharedgnms,],msvoom$weights[sharedgnms,])%>%
		set_colnames(alldesign$sample_id)
	allvoom$design<-model.matrix(data=alldesign,~group*assay)%>%
		set_rownames(alldesign$sample_id)
	allvoom$targets<-bind_rows(countvoom$targets,msvoom$targets)
	allvoom
}
stopifnot(exists('allvoom'))
stopifnot(exists('alldesign'))
stopifnot(exists('sharedgnms'))
stopifnot(exists('countvoom'))

################################################################################
########## Plot the variances
################################################################################
plotfile <- here(paste0('plots/countmeanvariance.pdf'))
pdf(plotfile)
limma::plotSA(lmFit(countvoom[sharedgnms,]), 
	main="Final model: Mean-variance trend", 
	ylab = "Sqrt( standard deviation )")
dev.off()
message(normalizePath(plotfile))
#
plotfile <- here(paste0('plots/msmeanvariance.pdf'))
pdf(plotfile)
plotSA(msfit, 
	main="Final model: Mean-variance trend", 
	ylab = "Sqrt( standard deviation )")
dev.off()
message(normalizePath(plotfile))


################################################################################
########## Use proDA to fit a linear model describing RNA-protein differences
###############################################################################
#Use the voom object's expression levels to fit a proDA model
stopifnot(exists('allvoom'))
stopifnot(exists('alldesign'))
ms_mat <- ms_data%>%{set_rownames(as.matrix(ms_data[,-1]),ms_data$gene_name)}
ms_mat[!is.finite(ms_mat)] <- NA
if(!file.exists(here('data/proDAfitms.rds'))){
	proDAfitms <- proDA::proDA(
			data=allvoom$E, 
			design = ~ group+assay+group:assay, 
            col_data = alldesign,
            n_subsample=1e3,
            reference_level = 'MB'
	)
	saveRDS(proDAfitms,here('data/proDAfitms.rds'))
}else{
	proDAfitms<-readRDS(here('data/proDAfitms.rds'))
}
#Now use proDA to get log fold changes corresponding to our contrasts
proda_res<- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	#get the columns for our  contrast
	grp1cols = colnames(allvoom)%>%str_subset(paste0(sampgroup_pair[[1]],'_\\d'))
	grp2cols = colnames(allvoom)%>%str_subset(paste0(sampgroup_pair[[2]],'_\\d'))
	g1ms <- grp1cols%>%str_subset('MS')%>%tail(1)
	g2ms <- grp2cols%>%str_subset('MS')%>%tail(1)
	g1rnaseq <- grp1cols%>%str_subset('rnaseq')%>%tail(1)
	g2rnaseq <- grp2cols%>%str_subset('rnaseq')%>%tail(1)
	#create numeric vectors which specify the contrast between our two samples in the
	#linear space of proDA's model
	contrlist = list(
		dev=(design(proDAfitms)[g2ms,] - design(proDAfitms)[g1ms,]) - 
			(design(proDAfitms)[g2rnaseq,] - design(proDAfitms)[g1rnaseq,])		,
		ms=(design(proDAfitms)[g2ms,] - design(proDAfitms)[g1ms,]),
		rna=(design(proDAfitms)[g2rnaseq,] - design(proDAfitms)[g1rnaseq,])
	)
	#now iterate over these to get fold changes for dev, ms and rna
	map_df(.id='effect',contrlist,function(contr){	
	proDA::test_diff(proDAfitms,contr=contr)%>%
		mutate(sig = pval<0.05)%>%
		select(gene_name=name,logFC=diff,sig)%>%
		mutate(contrast=contrnm)
		})
})
proda_res%<>%bind_rows

################################################################################
########## Finally produce quadrant plots with proDA
################################################################################
# contrnm=names(sampgroup_pairs)[1]
# sampgroup_pair=sampgroup_pairs[[1]]
#iterate over the pairs of samples and their names, plotting and defining
#our groups
quadrantplots <- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	stopifnot(c()%in%colnames(proda_res))
	compdata = proda_res%>%
		filter(contrast==contrnm)%>%
		pivot_wider(names_from=c('effect'),values_from=c('logFC','sig'))
	#define the classes
	compdata<-compdata%>%
	mutate(
		sigposrna =  sig_rna & (logFC_rna>0),
		signegrna =  sig_rna & (logFC_rna<0),
		nosigrna = !sig_rna,
		nosigdelta = !sig_dev,
		delta_pos = sig_dev & (logFC_dev > 0),
		delta_neg = sig_dev & (logFC_dev < 0),
		set = case_when(
		(sigposrna)&(nosigdelta) ~ 'RNA-Up',
		(sigposrna)&(delta_neg)&(abs(logFC_rna)<abs(logFC_dev)) ~ 'Divergent Down',
		(sigposrna)&(delta_neg) ~ 'RNA-Up, delta_buff',
		(sigposrna)&(delta_pos) ~ 'RNA-Up, delta_syn',
		(signegrna)&(nosigdelta) ~ 'RNA-down',
		(signegrna)&(delta_pos)&(abs(logFC_rna)<abs(logFC_dev)) ~ 'Divergent Up',
		(signegrna)&(delta_pos) ~ 'RNA-down, delta_buff',
		(signegrna)&(delta_neg) ~ 'RNA-down, delta_syn',
		(nosigrna)&(delta_pos) ~ 'Protein Up',
		(nosigrna)&(delta_neg) ~ 'Protein Down',
		TRUE ~ 'No Sig Effect'
	))
	#create plot
	plot = compdata%>%
		filter(!is.na(sig_rna),!is.na(sig_dev))%>%
		ggplot(aes(logFC_rna,logFC_ms,color=set,alpha=set!='No Sig Effect'))+
		geom_point(size=I(1),alpha=I(0.5))+
		scale_x_continuous(limits=c(-3,3))+
		scale_y_continuous(limits=c(-3,3))+
		theme_bw()+
		ggtitle(contrnm)+
		guides(color = guide_legend(override.aes = list(size=5)))
	#now print the plot	
	plotfile <- here(paste0('plots/','quadrant_lfc_',contrnm,'.pdf'))
	pdf(plotfile,w=10,h=10)
	print(plot)
	dev.off()
	message(normalizePath(plotfile))
	#
	#export
	tablnm<-here(paste0('tables/',contrnm,'_lfc_sets.tsv'))
	compdata%>%write_tsv(tablnm)
})
