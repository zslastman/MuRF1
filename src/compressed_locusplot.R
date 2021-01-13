#function that takes all the granges in a region, and compresses them down so introns
#are all a fixed distance
library(Gviz)
library(rtracklayer)
library(tidyverse)
library(magrittr)
library(Rsamtools)
library(GenomicAlignments)

grpcols = c(MB='#000099',MT='#FFFF00',MT_Ctrl24='#009900',MT_Ctrl72='#00FF00',MT_A24='#990000',MT_A72='#FF0000')

base::source(here::here('src/Rprofile.R'))

options(ucscChromosomeNames=FALSE) 

#this is the function that shifts a bunch of ranges to minimize the space between them, down to a ceiling maxwidth
maxwidth=1
get_comp_shift<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%.[-1]
	rangegaps%>%width
	toshift<-rep(0,length(rangestoshift))
	i=3
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			toshift[rightofgap] <- toshift[rightofgap] + adj
		}
	}
	toshift %<>% {. + min(start(rangestoshift))-1}
	#rangestoshift<-GenomicRanges::shift(rangestoshift,-rangestoshift$toshift)
	toshift
}

# annotation%<>%GenomicRanges::shift(-annotation$toshift)
totbams = Sys.glob('pipeline/star/data/*/*.bam')
totbams%<>%setNames(.,basename(dirname(.)))


# if(!exists('libsizes'))

if(TRUE){
	# peptidemsdata = fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/MS_Data_New/Ages_Brain_PEP_summ.txt')%>%
	# rename('gene_name'=Gene.names)

	# peptidemsdata%>%colnames
	# normpepmat = peptidemsdata%>%select(matches('Intensity.*input'))%>%as.matrix%>%
		# proDA::median_normalization(.)

	# normpepmat = colnames(normpepmat)%>%str_extract('(?<=\\.).*?(?=_input)')%>%
		# str_replace('p','')%>%
		# {split(seq_along(.),.)}%>%
		# lapply(.,function(colgrp){
			# normpepmat[,colgrp]%>%rowMedians(na.rm=T)
		# })%>%
		# simplify2array
	# rownames(normpepmat) <- peptidemsdata$Sequence

	extrseq = function(gr)GenomicFeatures::extractTranscriptSeqs(gr,x=FaFile('pipeline/GRCm38.p5.genome.chr_scaff.fa'))

	libstg = c(totbams)%>%dirname%>%basename%>%str_replace('_\\d+$','')
	libstg%<>%setNames(c(totbams)%>%dirname%>%basename)
	iso_tx_countdata <- readRDS(here('data/iso_tx_countdata.rds'))
	libsizes = iso_tx_countdata$counts%>%colSums%>%divide_by(1e6)

	# proteinsequences = GenomicFeatures::extractTranscriptSeqs(cdsgrl,x=FaFile('pipeline/GRCm38.p5.genome.chr_scaff.fa'))
	# aaproteinsequences = translate(proteinsequences)
	# trlist = ids_nrgname%>%distinct(gene_name,transcript_id)%>%{split(.[[2]],.[[1]])}
	# iso_tx_countdata <- readRDS(here('data/iso_tx_countdata.rds'))

}
gtf = here('ext_data/gencode.vM12.annotation.gtf')%>%normalizePath(mustWork=T)

# if(!exists('geneanno'))geneanno = fread(str_interp('grep -ie ${igene} ${gtf}'))%>%
if(!exists('geneanno'))geneanno = fread(gtf)%>%
	as.matrix%>%apply(1,paste0,collapse='\t')%>%paste0(collapse='\n')%>%import(text=.,format='gtf')

igene='Flna'

COMPRESS=T
# for(igene in c('Satb2','Flna','Nes','Bcl11b','Tle4'))
trimids <- function(x) x%>%str_replace('\\.\\d+','')
genesofinterest = fread('ext_data/2000929_TranscriptomeCandList.txt')%>%
	set_colnames(c('UniprotID','transcript_id'))

tr2gnm<-mcols(geneanno)%>%as.data.frame%>%mutate(transcript_id=transcript_id%>%trimids)%>%select(transcript_id,gene_name)
genesofinterest%>%left_join(tr2gnm,by='transcript_id')
genesofinterest_manual<-fread('src/genesofinterest.txt')
stopifnot(genesofinterest_manual$UniprotIDs%>%str_split(',')%>%unlist%>%is_in(genesofinterest[[1]])%>%all)

genesofinterest_manual$gene_name%>%intersect(geneanno$gene_name%>%unique)
genesofinterest_manual$gene_name%>%setdiff(geneanno$gene_name%>%unique)

genes2plot = genesofinterest_manual$gene_name
# genes2plot = c(genes2plot,c('Myh2','Myh6'))
genes2plot = c(c('Myh2','Myh6'))

for(COMPRESS in c(TRUE,FALSE)){
# for(igene in 'Trim63'){
for(igene in genes2plot){
{

types = c('UTR','CDS','exon')
if(!COMPRESS) types = c(types,'gene')

geneanno$gene_name%>%str_subset(igene)

igeneanno <- geneanno%>%
	subset(gene_name==igene)%>%
	subset(type%in%types)
stopifnot(length(igeneanno)>0)
if(COMPRESS){
	igeneanno$toshift <- get_comp_shift(igeneanno,maxwidth=10)
}else{
	igeneanno$toshift <- 0
}
# igeneanno$transcript_id%<>%trimids

#remove the exons for coding genes
codingtrs = igeneanno%>%subset(type=='CDS')%>%.$transcript_id
igeneanno%<>%subset(!((type=='exon')&(transcript_id%in%codingtrs)))


exontrack =
  igeneanno%>%
  GenomicRanges::shift(.,-.$toshift)%>%
  # {.$transcript_id %<>% str_replace('\\.[0-9]+$','');.}%>%
  subset(type%in%setdiff(types,'gene'))%>%
  .[,c('type','transcript_id')]%>%
  {.$feature<-as.character(.$type);.}%>%
  {.$transcript<-as.character(.$transcript_id);.}%>%
  Gviz::GeneRegionTrack(.,transcriptAnnotation='transcript',geneSymbol=FALSE,showId=TRUE,thinBoxFeature=c("exon","UTR"))
  
library(Gviz)
}


get_cov_track <- function(igeneanno,bams=ribobams,offsets_=offsets){
	#get the reads
	# igeneanno = igeneanno%>%subset(gene_name=igene)
	st = as.character(strand(igeneanno[1]))
	ichr = as.character(seqnames(igeneanno))[[1]]
	# psites = lapply(ribobams,function(bam)get_genomic_psites(bam,igeneanno,offsets_,comps=NULL))
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=200,which=igeneanno)
  	psites <- lapply(bams,function(bam)readGAlignments(bam,param=riboparam))

	psites%<>%lapply(coverage)
	normcov = map2(psites,libsizes[names(psites)],~ (.x / .y)*1e6)
	tpnormcov = normcov%>%{split(.,libstg[names(.)])} %>%lapply(function(trackpair) purrr::reduce(.f=`+`,trackpair))
	tpnormcov = tpnormcov%>%map(GRanges)%>%GRangesList%>%unlist%>%subset(score!=0)
	bp1tpnormcov=tpnormcov@ranges%>%as('IntegerList')%>%setNames(.,seq_along(.))%>%stack%>%as.data.frame
	bp1tpnormcov$grp = names(tpnormcov)[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov$score = tpnormcov$score[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov = bp1tpnormcov%>%select(start=value,grp,score)%>%
		spread(grp,score)%>%
		mutate(seqnames=ichr,end=start)%>%
		GRanges(.)
	

	ov = findOverlaps(bp1tpnormcov,igeneanno)
	ov = ov[!duplicated(queryHits(ov)),]
	bp1tpnormcov = GenomicRanges::shift(bp1tpnormcov[ov@from],-igeneanno$toshift[ov@to])
	strand(bp1tpnormcov) = st
	# bp1tpnormcov
	stopifnot(all(libstg%in%colnames(mcols(bp1tpnormcov))))
	# for(tp in tps) mcols(bp1tpnormcov)[[tp]]%<>%{log2(.+(.5*min(.,na.rm=T)))}
	DataTrack(bp1tpnormcov,type='hist',col = grpcols,groups=names(grpcols))
}
# get_cov_track(igeneanno)

 
# peptidemsdata_=peptidemsdata;cdsgrl_=cdsgrl;aaproteinsequences_ = aaproteinsequences
get_gene_pepgr = function(igene,peptidemsdata_=peptidemsdata,cdsgrl_=cdsgrl,aaproteinsequences_ = aaproteinsequences){
	gpep_df = peptidemsdata_%>%filter(gene_name==igene)
	peptides=gpep_df$Sequence
	# peptides = 'AAPAETDQR'
	ig_seqs = aaproteinsequences_[trlist[[igene]]]
	#
	peptidegr = peptides%>%setNames(.,.)%>%map_df(.id='peptide',function(peptide){
		vmatchPattern(peptide,ig_seqs)%>%
		as.data.frame%>%
		transmute(start,end,width,seqnames=names(ig_seqs)[group])
	})%>%GRanges

	gpeptidegr = peptidegr%>%{
		x=.
		x = resize(x,1,'start',ignore.strand=TRUE)
		end(x)  = ((start(x)-1)*3)+1
		start(x)  = end(x)
		# end(x)  = start(x)
		out = mapFromTranscripts(x,cdsgrl_)
		out$transcript_id = names(cdsgrl_)[out$transcriptsHits]
		out$peptide = x$peptide[out$xHits]
		out
		}
	# gpep	
	gpeptidegr$xHits=NULL
	gpeptidegr$transcriptsHits=NULL
	mcols(gpeptidegr)%<>%cbind(normpepmat[match(gpeptidegr$peptide,rownames(normpepmat)),])
	gpeptidegr
}
# get_gene_pepgr('Satb2')

get_peptide_track<-function(igeneanno){
	igeneanno%<>%subset(type!='gene')
	igene = igeneanno$gene_name[1]
	pepgr = get_gene_pepgr(igene)
	pepgr = pepgr[!duplicated(pepgr$peptide),]
	# igeneanno = annotation%>%subset(gene_name=igene)
	ov = findOverlaps(pepgr,igeneanno)
	ov = ov[pepgr$transcript_id[queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
	ov = ov[!duplicated(queryHits(ov)),]
	pepgr = GenomicRanges::shift(pepgr[ov@from],-igeneanno$toshift[ov@to])
	# start(pepgr)%>%txtplot
	pepgr$transcript_id=NULL
	pepgr = unique(pepgr)
	pepgr = sort(pepgr)
	# pepgr[60:70]
	# pepgr = pepgr[70:80]
	DataTrack(pepgr[,tps],cols = grpcols,groups=tps)	
	# pepgr
}

# motifpattern='TGTANATA'

# getmotiftrack = function(motifpattern,igeneanno){
# 		igeneanno%<>%subset(type!='gene')
# spligeneanno = igeneanno%>%split(.,.$transcript_id)
# trs = names(spligeneanno)
# igeneannoseq = spligeneanno%>%extrseq

# motifgr = vmatchPattern(motifpattern,fixed=F,igeneannoseq)%>%
# 	as.data.frame%>%
# 	transmute(start,end,width,seqnames=names(igeneannoseq)[group])%>%
# 	GRanges%>%
# 	mapFromTranscripts(spligeneanno)
# ov = findOverlaps(motifgr,igeneanno)
# ov = ov[trs[motifgr$transcriptsHits][queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
# ov = ov[!duplicated(queryHits(ov)),]
# motifgr = GenomicRanges::shift(motifgr[ov@from],-igeneanno$toshift[ov@to])
# motifgr = unique(motifgr)
# AnnotationTrack(motifgr)
# }


# getstart_track = function(igeneanno){
# 	spligeneanno = igeneanno%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%.[str_order_grl(.)]
# 	spligeneanno = spligeneanno%>%resize_grl(1,'start')%>%unlist
# 	spligeneanno %<>% GenomicRanges::shift(-spligeneanno$toshift)
# 	AnnotationTrack(spligeneanno,background.panel = "grey")
# }

# get_peptide_track('Satb2',annotation)

# ranges('Satb2',annotation)


#run plotting function
{
nametrack = function(x,tnm) x%>%{.@name=tnm;.}
#now plot
options(ucscChromosomeNames=FALSE)
compressstr = if(COMPRESS) '' else '_NONcompressed_'
plotfile<- here(paste0('plots/locusplot_',compressstr,igene,'.pdf'))
pdf(plotfile,w=24,h=12)
plotTracks(list(
	get_cov_track(igeneanno,totbams)%>%nametrack('RNA-Seq (RPM)'),
	# get_peptide_track(igeneanno)%>%nametrack('log2(Normalized Intensity)'),
 	exontrack%>%nametrack('Transcripts')
 	
),col = grpcols)
dev.off()
message(normalizePath(plotfile))
}

trid2gid = 'data/trid2gid.txt'%>%fread

iddf = 'data/gid_2_gname.txt'%>%
	fread%>%
	set_colnames(c('gnm','g_id'))%>%
	left_join(trid2gid)

dsets = iso_tx_countdata$abundance%>%colnames
colgrps = dsets%>%str_extract('.*?(?=_\\d+)')
colgrp = colgrps[[1]]

ucolgrps = unique(colgrps)%>%setNames(.,.)
relativeusedf = map_df(.id='dset',ucolgrps,function(colgrp){
	sab = iso_tx_countdata$abundance[,colgrps==colgrp]%>%rowSums
	sab%>%
		enframe('tr_id','abundance')%>%
		left_join(iddf%>%select(tr_id,gnm))%>%
		group_by(gnm)%>%mutate(abundance=abundance/sum(abundance))
})


ign_usedf = relativeusedf%>%
  safe_filter(gnm==igene)%>%
  # separate(dset,c('time','assay'))
  identity

#now plot
dir.create(here('plots/dtuplots'))

plotfile<- here(paste0('plots/dtuplots/',igene,'_DTUplot','.pdf'))
cairo_pdf(plotfile)
print(ign_usedf%>%
  ggplot(data=.,aes(x=dset,fill=tr_id,y=abundance))+
  geom_bar(stat='identity')+
  # scale_fill_discrete(name='colorname',colorvals)+
  # facet_grid(assay~.)+
  scale_x_discrete(paste0('sample group'))+
  scale_y_continuous(paste0('relative usage'))+
  ggtitle(paste0('Transcript Usage ',igene))+
  theme_bw())
dev.off()
message(normalizePath(plotfile))



}
}

