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
#
################################################################################
########## t tests
################################################################################


sampgroup_pair<-sampgroup_pairs[[1]]
ttestdfs <- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	# contrnm<-names(sampgroup_pairs)[[1]]
	grp1 = sampgroup_pair[[1]]
	grp2 = sampgroup_pair[[2]]
	#
	grp1cols = colnames(ms_data)%>%str_subset(paste0(grp1,'_\\d'))
	grp2cols = colnames(ms_data)%>%str_subset(paste0(grp2,'_\\d'))
	#
	ttestlist <- function(list){
		require(broom)
		tidy(t.test(list[[1]],list[[2]]))
	}
	ttestlist<-possibly(ttestlist,NULL)
	#
	msttests<-ms_data%>%
		select(gene_name,one_of(c(grp1cols,grp2cols)))%>%
		pivot_longer(-gene_name,names_to='dataset',values_to='value')%>%
		mutate(grp = dataset%>%str_replace('_\\d+$',''))%>%
		# filter(gene_name==testgene)
		group_by(gene_name)%>%
		summarise(ttest = 
			list(value[dataset%in%grp2cols],value[dataset%in%grp1cols])%>%
			ttestlist%>%list
		)
	msttests%>%sample_n(10)
	#
	msttests%<>%filter(map_lgl(ttest,Negate(is.null)))
	#
	msttests%<>%unnest(ttest)%>%
		mutate(p.value=p.adjust(p.value))
	msttests%<>%mutate(contrast = contrnm)
	msttests
})
ttestdfs%<>%bind_rows

grp1='MB'
grp1cols = colnames(ms_data)%>%str_subset(paste0(grp1,'_\\d'))
plotfile <- here(paste0('plots/','zeros_vs_mean','.pdf'))
pdf(plotfile)
qplot(
	x=ms_data[,grp1cols]%>%apply(1,.%>%keep(is.finite)%>%median),
	y=ms_data[,grp1cols]%>%apply(2,Negate(is.finite))%>%rowSums,
	geom='bin2d'
)+geom_smooth()+geom_jitter()
dev.off()
message(normalizePath(plotfile))





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
#
# subsetsamps <- c(grp1cols,grp2cols)%>%paste0('MS_',.)
#
subsetvoom<-function(subsetsamps, allvoom){
	subvoom<-new('EList')
	subvoom$E<-allvoom$E[,subsetsamps]
	subvoom$weights<-allvoom$weights[,subsetsamps]
	subvoom$design<-allvoom$design[subsetsamps,]
	subvoom$targets<-allvoom$targets[subsetsamps,,drop=F]
	descols<-subvoom$design%>%colSums%>%.[.!=0]%>%names
	subvoom$design%<>%.[,descols]
	subvoom	
}
#
}


################################################################################
########## Get limma based mass spec estimates
################################################################################
{
sampgroup_pair <- sampgroup_pairs[[1]]
#
#limma ms contrasts
mslimmacontrasts <- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	msvoom <- msvoom
	grp1 = sampgroup_pair[[1]]
	grp2 = sampgroup_pair[[2]]
	grp1cols = colnames(msvoom$E)%>%str_subset(paste0(grp1,'_\\d'))
	grp2cols = colnames(msvoom$E)%>%str_subset(paste0(grp2,'_\\d'))
	contrcols <- c(grp1cols,grp2cols)
	contrvoom <- subsetvoom(contrcols,msvoom)
	contrvoom$design<-msdesign%>%
		filter(sample_id%in%contrcols)%>%
		arrange(group==grp2)%>%
		mutate(group = as_factor(as.character(group)))%>%
		model.matrix(data=.,~group)
	contrastebayes <- contrvoom%>%lmFit%>%eBayes
	contrcoef = contrvoom$design%>%colnames%>%str_subset(grp2)
	message(contrcoef)
	topTable(contrastebayes,n=Inf,coef=contrcoef)%>%
		as.data.frame%>%
		mutate(gene_name=rownames(.))%>%
		set_rownames(NULL)%>%
		# filter(gene_name==testgene)
		identity
})
#
mslimmacontrasts%<>%bind_rows(.id='contrast')
# 
# #check to compare ttests and limma
# mslimmacontrasts%>%inner_join(ttestdfs)%>%
# 	filter(contrast=='MB_MT')%>%
# 	# {lm(data=.,estimate~logFC)}%>%
# 	# {txtplot(.$estimate,.$logFC)}%>%
# 	{cor.test(.$estimate,.$logFC)}

# ttestdfs%>%filter(contrast=='MB_MT')%>%filter(p.value<0.05)%>%arrange(estimate)


# # mslimmacontrasts%>%filter(contrast=='MB_MT')%>%
# # 	filter(gene_name==testgene)
# # ttestdfs%>%filter(contrast=='MB_MT')%>%
# # 	filter(gene_name==testgene)

# ms_data%>%
# 	select(gene_name,one_of(c(grp1cols,grp2cols)))%>%
# 	pivot_longer(-gene_name,names_to='dataset',values_to='value')%>%
# 	filter(gene_name==testgene)%>%
# 	group_by(dataset%in%grp2cols)%>%
# 	summarise(mean(value))

# ((msvoom$E[testgene,grp2cols]%>%mean)- (msvoom$E[testgene,grp1cols]%>%mean))
# ((allvoom$E[testgene,grp2cols]%>%mean)- (allvoom$E[testgene,grp1cols]%>%mean))
# ((allvoom$E[testgene,grp2cols]%>%mean)- (allvoom$E[testgene,grp1cols]%>%mean))
# alldesign
# allvoom$E%>%colnames
# (ms_data%>%filter(gene_name==testgene)%>%.[,grp2cols]%>%unlist%>%mean)- 
# 	(ms_data%>%filter(gene_name==testgene)%>%.[,grp1cols]%>%unlist%>%mean)
}

################################################################################
########## Now get the log fold changes for the quadrant plots
################################################################################
{
contrnm='MT_Ctrl24_MT_A24'
sampgroup_pair=sampgroup_pairs[[contrnm]]
alllimmacontrasts <- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	grp1 = sampgroup_pair[[1]]
	grp2 = sampgroup_pair[[2]]
	contrcols = colnames(allvoom$E)%>%
		str_subset(paste0(grp1,'_\\d','|',grp2,'_\\d'))
	contrcols%<>%.[order(str_detect(.,'_A'))]
	contrcols%<>%.[order(!str_detect(.,'rnaseq'))]
	contrvoom <- subsetvoom(contrcols,allvoom)
	contrvoom$E%>%colnames
	# (contrvoom$E[testgene,grp2cols]%>%mean)-(contrvoom$E[testgene,grp1cols]%>%mean)
	# contrvoom$E[1,]<-c(
	# 	1,1,1,1,1*2,1*2,1*2,1*2,
	# 	2,2,2,2,2,2,2,2,2,2
	# )
	# (contrastebayes$design%*%contrastebayes$coef[testgene,])[grp2cols,]
	# (contrastebayes$design%*%contrastebayes$coef[testgene,])[grp1cols,]
	contrvoom$design<-alldesign%>%
		filter(sample_id%in%contrcols)%>%
		arrange(group==grp2)%>%
		mutate(group = as_factor(as.character(group)))%>%
		{set_rownames(model.matrix(data=.,~group*assay),.$sample_id)}%>%
		.[contrvoom$E%>%colnames,]
	contrastebayes <- contrvoom%>%lmFit%>%eBayes
	coefs = contrvoom$design%>%colnames%>%str_subset(neg=T,'Intercept')%>%
		setNames(.,.)
	contrastebayes$design
	contrvoom$E%>%colnames
	coefsout = map_df(.id='coef',coefs,function(coef){
		topTable(contrastebayes,n=Inf,coef=coef)%>%
			as.data.frame%>%
			mutate(gene_name=rownames(.))%>%
			set_rownames(NULL)
	})%>%mutate(contrast=contrnm)
	coefsout
})
alllimmacontrasts%<>%bind_rows
}

# #
# alllimmacontrasts%>%bind_rows%>%sample_n(10)
# #
# alllimmacontrasts%>%bind_rows%>%
# 	filter(contrast==testcontr)%>%
# 	group_by(coef)%>%slice(1)%>%
# 	as.data.frame
# #
# txnchange<-alllimmacontrasts%>%
# 	filter(contrast==testcontr)%>%
# 	# filter(gene_name==gene_name[1])%>%
# 	filter(coef%>%str_detect(':'))
# inteff<-alllimmacontrasts%>%
# 	filter(contrast==testcontr)%>%
# 	# filter(gene_name==gene_name[1])%>%
# 	filter(coef%>%str_detect('^[^:]+$'))%>%
# 	filter(coef%>%str_detect(neg=T,'assay'))
# #
# left_join(txnchange,inteff,by='gene_name')%>%{txtplot(.$logFC.x,.$logFC.y)}
# left_join(txnchange,inteff,by='gene_name')%>%{txtplot(.$logFC.x,.$logFC.x+.$logFC.y)}

# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%mutate(contrast=contrnm)%>%
# 	left_join(.,inteff,by='gene_name')%>%{txtplot(.$logFC.x,.$logFC.y)}

# testcontr<-testcontr

# left_join(txnchange,
# 	mslimmacontrasts%>%filter(contrast==testcontr),by=c('gene_name'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}


# left_join(txnchange,
# 	mslimmacontrasts%>%filter(contrast==testcontr),by=c('gene_name'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}

# #okay so estimation of the mass spec fold changes works
# mslfc%>%enframe('gene_name','logFC')%>%mutate(contrast=contrnm)%>%
# 	inner_join(mslimmacontrasts,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}

# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%mutate(contrast=contrnm)%>%
# 	inner_join(mslimmacontrasts,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}

# #
# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%
# 	mutate(contrast=contrnm)%>%
# 	inner_join(ttestdfs,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC,.$estimate)}%>%
# 	{cor.test(.$logFC,.$estimate)}

# alllimmacontrasts%>%
# 	filter(contrast==testcontr)%>%
# 	filter(coef%>%str_detect('group[^:]+$'))%>%
# 	inner_join(mslimmacontrasts,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}



plotfile <- here(paste0('plots/','deltaMS_quadrants_',contrnm,'.pdf'))
pdf(plotfile)
p=alllimmacontrasts%>%
	filter(contrast==contrnm)%>%
	# filter(gene_name==testgene)%>%
	mutate(sig = adj.P.Val<0.05)%>%
	filter(coef!='assayMS')%>%
	mutate(txncoef=coef%>%str_detect('group[^:]+$'))%>%
	mutate(coef=ifelse(!txncoef,'deltaFC_MS','Shared'))%>%
	select(coef,sig,logFC,gene_name)%>%
	pivot_wider(names_from='coef',values_from=c('sig','logFC'))%>%
	ggplot(data=.,aes(x=logFC_Shared,y=logFC_deltaFC_MS))+
	geom_point()+
	theme_bw()
print(p)
dev.off()
message(normalizePath(plotfile))

################################################################################
########## Plot the variances
################################################################################
plotfile <- here(paste0('plots/countmeanvariance.pdf'))
pdf(plotfile)
plotSA(countfit[sharedgnms,], 
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
########## Perhaps just use proDA
################################################################################

ms_mat <- ms_data%>%{set_rownames(as.matrix(ms_data[,-1]),ms_data$gene_name)}
ms_mat[!is.finite(ms_mat)] <- NA
if(!file.exists(here('data/proDAfitms.rds'))){
	proDAfitms <- proDA::proDA(
			data=allvoom$E, 
			design = ~ group+assay+group:assay, 
			# ms_mat[,msdesign$sample_id], 
			# design = ~ group+assay+group:assay, 
			# design = ~ group, 
            col_data = alldesign,
            # col_data = msdesign,
            n_subsample=1e3,
            # reference_level = "E13"
            reference_level = 'MB'
	)
	saveRDS(proDAfitms,here('data/proDAfitms.rds'))
}else{
	proDAfitms<-readRDS(here('data/proDAfitms.rds'))
}


proda_res<- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
	#
	grp1cols = colnames(allvoom)%>%str_subset(paste0(sampgroup_pair[[1]],'_\\d'))
	grp2cols = colnames(allvoom)%>%str_subset(paste0(sampgroup_pair[[2]],'_\\d'))
	g1ms <- grp1cols%>%str_subset('MS')%>%tail(1)
	g2ms <- grp2cols%>%str_subset('MS')%>%tail(1)
	g1rnaseq <- grp1cols%>%str_subset('rnaseq')%>%tail(1)
	g2rnaseq <- grp2cols%>%str_subset('rnaseq')%>%tail(1)
	#
	contrlist = list(
		dev=(design(proDAfitms)[g2ms,] - design(proDAfitms)[g1ms,]) - 
			(design(proDAfitms)[g2rnaseq,] - design(proDAfitms)[g1rnaseq,])		,
		ms=(design(proDAfitms)[g2ms,] - design(proDAfitms)[g1ms,]),
		rna=(design(proDAfitms)[g2rnaseq,] - design(proDAfitms)[g1rnaseq,])
	)
	map_df(.id='effect',contrlist,function(contr){	
	proDA::test_diff(proDAfitms,contr=contr)%>%
		mutate(sig = pval<0.05)%>%
		select(gene_name=name,logFC=diff,sig)%>%
		mutate(contrast=contrnm)
		})
})
proda_res%<>%bind_rows

# proda_res%>%group_by(gene_name,contrast,effect)%>%filter(n()>1)%>%group_slice(1)
# proda_res%>%filter(gene_name=='Gnai3')%>%arrange(contrast,effect)
# proda_res%>%filter(contrast=='MT_Ctrl24_MT_Ctrl72',gene_name=='Gnai3')


# contvect <- result_names(proDAfitms)%>%is_in(c('groupMT','`groupMT:assayMS`'))%>%
# 	as.numeric
# prodacontrdf<-proDA::test_diff(proDAfitms,contr=contvect)%>%
# 	select(gene_name=name,logFC=diff)%>%mutate(contrast=contrnm)
# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%
# 	mutate(contrast=contrnm)%>%
# 	inner_join(prodacontrdf,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}


dds <- readRDS('data/dds.rds')
design(dds) <- as.formula('~group')
dds <- DESeq(dds)
resnames <- c("Intercept", "group_MT_vs_MB", "group_MT_A24_vs_MB", "group_MT_A72_vs_MB", 
"group_MT_Ctrl24_vs_MB", "group_MT_Ctrl72_vs_MB")
stopifnot(resultsNames(dds)==resnames)


rnaseqlfc%>%
	inner_join(proda_res%>%filter(effect=='rna',contrast==contrnm),by=c('gnm'='gene_name'))%>%
	filter(log2FoldChange>-10)%>%
	{txtplot(.$logFC,.$log2FoldChange)}

rnaseqlfc%>%filter(log2FoldChange>3)%>%head

highgenes = rownames(tx_countdata$counts)[tx_countdata$counts%>%rowMin%>%`>`(10)]
tgene = rnaseqlfc%>%
	filter(gnm %in% highgenes)%>%
	filter(log2FoldChange>3)%>%.$gnm%>%.[1]

gnm2gid = setNames(names(gid2gnm),gid2gnm)

allvoom$E[tgene,]%>%.[names(.)%>%str_detect('MB|MT_\\d')]%>%
	enframe%>%separate(name,c('assay','cell','rep'))%>%
	group_by(assay,cell)%>%summarise_at(vars(value),mean)

normcountmat <- tx_countdata$counts%>%sweep(2,F='/',STATS=colSums(.))
normcountmat[tgene,]%>%.[names(.)%>%str_detect('MB|MT_\\d')]%>%
	enframe%>%separate(name,c('assay','cell','rep'))%>%
	group_by(assay,cell)%>%summarise_at(vars(value),mean)


################################################################################
########## Finally produce quadrant plots with proDA
################################################################################
contrnm=names(sampgroup_pairs)[1]
sampgroup_pair=sampgroup_pairs[[1]]
quadrantplots <- imap(sampgroup_pairs,function(sampgroup_pair,contrnm){
    	# samp1=sampgroup_pair[[1]]
     #  samp2=sampgroup_pair[[2]]
     #  msdatacols1 = ms_data%>%colnames%>%
     #    mscolnames_replace%>%
     #    str_detect(paste0(samp1,'_\\d'))
     #  msdatacols2 = ms_data%>%colnames%>%
     #    mscolnames_replace%>%
     #    str_detect(paste0(samp2,'_\\d'))
     #  stopifnot(any(msdatacols1))
     #  stopifnot(any(msdatacols2))
     #  msdata1 = ms_data[,msdatacols1]%>%rowMeans(na.rm=T)%>%
     #    setNames(rownames(ms_data))
     #  msdata2 = ms_data[,msdatacols2]%>%rowMeans(na.rm=T)%>%
     #    setNames(rownames(ms_data))
     #  mslfc = ((msdata2) - (msdata1))%>%
     #  	setNames(ms_data$gene_name)
     #  # mslfc%>%keep(is.finite)%>%txtdensity
     #  #

     #  rnaseqlfc = results(dds,c('group',samp2,samp1))
     #  #
     #  rnaseqlfc= rnaseqlfc%>%as.data.frame%>%
     #  	rownames_to_column('g_id')%>%
     #  	# select(-g)
     #    left_join(gid2gnm%>%enframe('g_id','gnm'),by='g_id')%>%
     #    # filter(padj<0.05)%>%
     #    mutate(rnaseq_sig = padj<0.05)%>%
     #    select(gnm,log2FoldChange,rnaseq_sig)
     #  # repgnms = rnaseqlfc$gnm%>%table%>%keep(~.>1)%>%names
     #  #
     #  compdata = rnaseqlfc%>%
     #  	mutate(contrast=contrnm)%>%
     #    inner_join(enframe(mslfc,'gnm','ms_l2fc'),by='gnm')%>%
     #    inner_join(proda_res,by=c('gnm'='gene_name','contrast'))
     #  compdata%>%colnames
      compdata = proda_res%>%
      	filter(contrast==contrnm)%>%
      	pivot_wider(names_from=c('effect'),values_from=c('logFC','sig'))
      #
      plotfile <- here(paste0('plots/','quadrant_lfc_',contrnm,'.pdf'))
      pdf(plotfile,w=10,h=10)
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
      plot = compdata%>%
      	filter(!is.na(sig_rna),!is.na(sig_dev))%>%
      	# filter(set%in%c('RNA-Up','RNA-down','RNA-Up, delta_buff','RNA-Up, delta_syn'))%>%
      	# mutate(agr=sign(log2FoldChange)==sign(logFC))%>%
      	ggplot(aes(logFC_rna,logFC_ms,color=set,alpha=set!='No Sig Effect'))+
      	geom_point(size=I(1),alpha=I(0.5))+
      	scale_x_continuous(limits=c(-3,3))+
      	scale_y_continuous(limits=c(-3,3))+
      	#facet_grid()+
      	# geom_text(data=tibble(labl=),
      		# aes(label=labl,x=Inf,y=Inf),hjust=1,vjust=1)+
      	theme_bw()+
      	ggtitle(contrnm)+
				guides(color = guide_legend(override.aes = list(size=5)))
			print(plot)
		dev.off()
      message(normalizePath(plotfile))
      #
      tablnm<-here(paste0('tables/',contrnm,'_lfc_sets.tsv'))
      compdata%>%
      	write_tsv(tablnm)
})

# RNAseq level changes.
# To compare changes at the transcriptional and protein level, we matched selected
# the Majority Protein ID class with the highest median signal, for each gene name
# We used limma ___ to remove a single batch effect for the data, and found that
# removing further batch effects did not improve the TPM-BAC correlation, nor did
# the use of MAD-normalisation techniques.

#To produce the quadrant plots in Figure ___, we made use of a dropout aware 
# linear-model, (proDA reference), in which log protein and variance stabilized
# (via limma's voom function) transcript expression were expresssed using common
# condition dependent effect (rna) and a protein specific 'MSdelta' effect. Based on the
# significance of this MSdelta effect, and the original DESeq2 results, We then 
# classified genes as follows:

# RNAseq significant positive, MSdelta not significant ~ 'RNA-down',
# RNAseq significant positive, MSdelta significant, positive ~ 'RNA-down, delta_syn',
# RNAseq significant positive, MSdelta significant, negative ~ 'RNA-down, delta_buff',
# RNAseq significant positive, MSdelta significant, negative, abs(MSdelta) > abs(RNA_l2fc) ~ 'Divergent Negative'
# RNAseq significant negative, MSdelta not significant ~ 'RNA-down',
# RNAseq significant negative, MSdelta significant, positive ~ 'RNA-down, delta_buff',
# RNAseq significant negative, MSdelta significant, negative ~ 'RNA-down, delta_syn',
# RNAseq significant negative, MSdelta significant, positive, abs(MSdelta) > abs(RNA_l2fc) ~ 'Divergent Positive'
# RNAseq no significant change, MSdelta significant, positive ~ 'Protein Up',
# RNAseq no significant change, MSdelta significant, negative ~ 'Protein Down',
# Otherwise:  'No Sig Effect'


# ################################################################################
# ########## 
# ################################################################################


# ###Try justlinear models?


# lmdfs <- imap(sampgroup_pairs[1],function(sampgroup_pair,contrnm){
# 	# contrnm<-names(sampgroup_pairs)[[1]]
# 	grp1 = sampgroup_pair[[1]]
# 	grp2 = sampgroup_pair[[2]]
# 	#
# 	grp1cols = colnames(allvoom)%>%str_subset(paste0(grp1,'_\\d'))
# 	grp2cols = colnames(allvoom)%>%str_subset(paste0(grp2,'_\\d'))
# 	#
# 	msttests<-allvoom$E%>%as.data.frame%>%rownames_to_column('gene_name')%>%
# 		select(gene_name,one_of(c(grp1cols,grp2cols)))%>%
# 		pivot_longer(-gene_name,names_to='dataset',values_to='value')%>%
# 		mutate(grp = dataset%>%str_replace('_\\d+$','')%>%str_replace('rnaseq_|MS_',''))%>%
# 		mutate(assay = dataset%>%str_extract('rnaseq|MS')%>%as_factor)%>%
# 		# filter(gene_name==testgene)%>%
# 		group_by(gene_name)%>%
# 		nest%>%
# 		summarise(ttest = map(data,~tidy(lm(data=.x,value~grp*assay))))%>%
# 		unnest(ttest)
# 	#
# 	# msttests%<>%filter(map_lgl(ttest,Negate(is.null)))
# 	#
# 	msttests%<>%mutate(p.value=p.adjust(p.value))
# 	msttests%<>%mutate(contrast = contrnm)
# 	msttests
# })
# lmdfs%<>%bind_rows

# lmcontdf<-lmdfs%>%filter(contrast==contrnm)%>%
# 	filter(term%>%str_detect(neg=T,'MS|Intercept'))%>%
# 	select(gene_name,logFC=estimate,contrast)

# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%
# 	mutate(contrast=contrnm)%>%
# 	inner_join(lmcontdf,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}

# lmcontdf<-lmdfs%>%
# 	filter(contrast==contrnm)%>%
# 	filter(term%>%str_detect(neg=T,'Intercept|^assayMS$'))%>%
# 	group_by(gene_name,contrast)%>%summarise(estimate=sum(estimate))%>%
# 	select(gene_name,logFC=estimate,contrast)
# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%
# 	mutate(contrast=contrnm)%>%
# 	inner_join(lmcontdf,by=c('gene_name'`,'contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}

# lmcontdf<-lmdfs%>%filter(contrast==contrnm)%>%
# 	filter(term%>%str_detect('\\:'))%>%
# 	group_by(gene_name,contrast)%>%summarise(estimate=sum(estimate))%>%
# 	select(gene_name,logFC=estimate,contrast)

# rnaseqlfc%>%set_colnames(c('gene_name','logFC'))%>%
# 	mutate(contrast=contrnm)%>%
# 	inner_join(lmcontdf,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}


# #okay so proDA and hte t-tests are for all intents and purposes, identical
# #wrt to total ms change
# lmcontdf<-lmdfs%>%
# 	filter(contrast==contrnm)%>%
# 	filter(term%>%str_detect(neg=T,'Intercept|^assayMS$'))%>%
# 	group_by(gene_name,contrast)%>%summarise(estimate=sum(estimate))%>%
# 	select(gene_name,logFC=estimate,contrast)
# # contvect <- result_names(proDAfitms)%>%is_in(c('groupMT','`groupMT:assayMS`'))%>%
# 	# as.numeric
# prodacontrdf<-proDA::test_diff(proDAfitms,contr=contvect)%>%
# 	select(gene_name=name,logFC=diff)%>%mutate(contrast=contrnm)
# lmcontdf%>%set_colnames(c('gene_name','logFC','contrast'))%>%
# 	inner_join(prodacontrdf,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}
# #and also the interaciton term
# lmcontdf<-lmdfs%>%
# 	filter(contrast==contrnm)%>%
# 	filter(term%>%str_detect('\\:'))%>%
# 	group_by(gene_name,contrast)%>%summarise(estimate=sum(estimate))%>%
# 	select(gene_name,logFC=estimate,contrast)
# contvect <- result_names(proDAfitms)%>%is_in(c('`groupMT:assayMS`'))%>%
# 	as.numeric
# prodacontrdf<-proDA::test_diff(proDAfitms,contr=contvect)%>%
# 	select(gene_name=name,logFC=diff)%>%mutate(contrast=contrnm)
# lmcontdf%>%set_colnames(c('gene_name','logFC','contrast'))%>%
# 	inner_join(prodacontrdf,by=c('gene_name','contrast'))%T>%
# 	{txtplot(.$logFC.x,.$logFC.y)}%>%
# 	{cor.test(.$logFC.x,.$logFC.y)}


# #proda estimates much more for the interaction term
# proDA::test_diff(proDAfitms,contr='`groupMT:assayMS`')%>%filter(adj_pval<0.05)%>%nrow
# lmdfs%>%filter(term%>%str_detect('\\:'))%>%filter(p.value<0.05)%>%nrow

# #and the shared term...
# proDA::test_diff(proDAfitms,contr='groupMT')%>%filter(adj_pval<0.05)%>%nrow
# lmdfs%>%filter(term%>%str_detect('grpMT'))%>%filter(p.value<0.05)%>%nrow



# #verbose up tehre... but in general the output of all these models seems to be a little
# #less correlated thant he output of good old fashioned t-tests., even lms.
# #So probably, we use proDA to call significance...

# #