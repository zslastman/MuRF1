uniprotmap = fread('ext_data/silac/silac_out.txt')%>%set_colnames(c('UniprotID','tr_ids'))
unmappeduniprot = 'ext_data/silac/unmapped_uniprot.txt'%>%fread%>%set_colnames(c('UniprotID'))%>%tail(-1)
uniprotfaids = fread('grep -e \'>\' ../cortexomics/ext_data/uniprot.MOUSE.2014-10.fasta',header=F)%>%mutate(mousename=str_extract(V3,'\\w+'))
silacdata = fread('ext_data/silac/20110913_Jida_2x4, r2_Diff_and_atrophy_für_Gunnar.txt')
silacdata = readxl::read_xlsx('ext_data/2021Jan_Results_Proteome.xlsx',sheet=1)

inclusiontable<-function(a,b){
	u = union(a,b)
	a=unique(a)
	b=unique(b)
	table(u%in%a,u%in%b)
}
# BiocManager::install(c('UniProt.ws'))

suppressPackageStartupMessages({
	library(UniProt.ws)
})

up <- UniProt.ws(taxId=10090)

uniprotkbs = keys(up, "UNIPROTKB")
uniparcs = keys(up, "UNIPARC")
unigids = keys(up, "GENEID")

uniprottrs = AnnotationDbi::select(up,uniprotkbs,select='UNIPROTKB',columns='ENSEMBL_TRANSCRIPT')

 # silacdata%>%colnames%>%head(10)
 # [1] "GeneNames"                   "ProteinID"
 # [3] "Uniprot"                     "Uniprot.Name"
 # [5] "Peptide.Sequence"            "MW.kDa"
 # [7] "Sequence.Length"             "Ratio.M.L.Normalized.Eluate"
 # [9] "Ratio.H.M.Normalized.Eluate" "Ratio.M.L.Normalized.Murf1"
silacdata%>%.$Uniprot.Name%>%is_in(egs)

unmappedsilacdata= silacdata%>%filter(
	Uniprot%>%str_split(';')%>%map(is_in,unmappeduniprot$UniprotID)%>%map_lgl(all)
)
idmatchdf = unmappedsilacdata %>% group_by(GeneNames,Uniprot)%>%transmute(
	gnamematch = GeneNames%>%str_split(';')%>%map(is_in,gid2gname$gnm)%>%map_lgl(any),
	unipmatch = Uniprot%>%str_split(';')%>%map(is_in,uniprotkbs)%>%map_lgl(any)
	# pidmatch = ProteinID%>%str_split(';')%>%map(is_in,uniparcs)%>%map_lgl(all),
)
idmatchdf%>%filter(!unipmatch)


#most gene names appear only once, a small percentage twice
silacdata%>%select(Uniprot,GeneNames)%>%mutate(GeneName = str_split(GeneNames,';'))%>%unnest(GeneName)%>%
		mutate(hasmatch = GeneName %in%gid2gname$gnm)%>%
		group_by(GeneNames)%>%mutate(n_match= sum(hasmatch))%>%
		filter(n_match==1)%>%
		tally


silacdata%>%select(Uniprot,GeneNames)%>%mutate(GeneNames = str_split(GeneNames,';'))%>%unnest(GeneNames)%>%
	.$GeneNames%>%table%>%table


silacdata%>%select(Uniprot,GeneNames)%>%mutate(Uniprot = str_split(Uniprot,';'))%>%unnest(Uniprot)%>%
	mutate(hasmatch = Unprot %in%egs)%>%
	group_by(GeneNames)%>%summarise(n_match=any(hasmatch))%>%
	.$hasmatch%>%table
	

silacdata$Uniprot%>%enframe%>%mutate(value = str_split(value,';'))%>%unnest(value)%>%
	.$value%>%setdiff(egs)

unmappeduniprot[[1]]%>%setdiff(egs)

uniprottrs = AnnotationDbi::select(up,
	unmappeduniprot$UniprotID,
	'ENSEMBL_TRANSCRIPT',
	'PDB'
	)

uniprottrs

unmappeduniprot[[1]]%>%is_in(egs)


#alsoget uniprotID-ensembl_peptide links from biomart
mousemart<- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
getBM<-mymemoise(getBM)
bm <- getBM(filters='uniprotswissprot',
	attributes=c('uniprotswissprot','ensembl_peptide_id'),values=unmappeduniprot[[1]],mart=mousemart)

#load the swissprot data from this gencode version
gc_protiddf<-fread(here('../cortexomics/annotation//gencode.vM12.metadata.SwissProt'),header=F)%>%set_colnames(c('transcript_id','uniprotkb_id','swissprot_id'))

#
#so all but 12 of those rows are in the uniptrobkb from ws objects
#800 of the unmappeduniprot ids are not even in the uniprotkbs object


# BiocManager::install('hiReadsProcessor')
library(hiReadsProcessor)

#blat -t=protein -q=protein ext_data/silac/unmapped_uniprot.fa   ../cortexomics/ext_data/gencode.vM12.pc_translations.fa  match.psl
unmappedpsl=read.psl('match.psl')
pslcols = 'matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts'%>%str_split(',')%>%.[[1]]
psldf = fread('match.psl',skip=3)%>%set_colnames(pslcols)
library(txtplot)
blatmatch = psldf%>%group_by(tName)%>%
	filter(matches==max(matches))%>%
	mutate(matchfrac = matches/tSize)%>%
	filter(matchfrac>.9)%>%
	select(UniprotID=tName,tr_ids=qName)%>%
	mutate(UniprotID=str_extract(UniprotID,'(?<=\\|)[^|]+'))%>%
	mutate(tr_ids=str_extract(tr_ids,'ENSMUST\\w+'))

allunimatch = bind_rows(blatmatch,uniprotmap)%>%group_by(UniprotID)%>%mutate(tr_id=str_split(tr_ids,','))%>%unnest(tr_id)
allunimatch%<>%select(-tr_ids)



################################################################################
########Easy way - just take inambig gene name matchesj
################################################################################
matchsilac = silacdata%>%filter(!str_detect(GeneNames,';'))%>%filter(GeneNames%>%is_in(gid2gname$gnm))

tx2gene <-  read_tsv(here('pipeline/gid_2_trid.txt.gz'))%>%set_colnames(c('g_id','tr_id'))

silacmatch = silacdata%>%select(Uniprot,GeneNames)%>%
	mutate(UniprotID = str_split(Uniprot,';'))%>%
	unnest(UniprotID)%>%
	left_join(allunimatch,by='UniprotID')%>%
	left_join(tx2gene,by='tr_id')%>%
	left_join(gid2gname,by='g_id')	

#how often does stuff map to more than one gene
silacmatch%>%group_by(Uniprot)%>%summarise(gnms = n_distinct(gnm))%>%.$gnms%>%table
silacmatchdf = silacmatch%>%group_by(Uniprot)%>%filter(n_distinct(gnm)==1)%>%distinct(gnm)

msilacdata = silacdata%>%inner_join(silacmatchdf)

msilacdata%>%colnames

msilacdata = msilacdata

###comapre MT to MB transition
{
mtmbsilac = msilacdata%>%select(gnm,lfc_MT_MB = Ratio.M.L.Normalized.Eluate,pval_MT_MB = Ratio.M.L.Significance.B..Eluate)

mrnachagne = resultslist$MT_vs_MB%>%select(gene_id,log2FoldChange,padj)%>%left_join(gid2gname,by=c(gene_id='g_id'))

mrna_silac_df = mtmbsilac%>%left_join(mrnachagne)

{
cortext = mrna_silac_df%>%filter(padj<0.05)%>%
# filter(gnm!='Hspd1')%>%
{cor.test(.$log2FoldChange,.$lfc_MT_MB)}%>%tidy%>%
	{str_interp('rho = ${round(.$estimate,3)},p = ${format(.$p.value,3,digits=2)}')}

#now plot
plotfile<- here(paste0('plots/','mrna_v_silac_fc_mt_mb','.pdf'))
pdf(plotfile)
print(
	mrna_silac_df%>%
	filter(padj<0.05)%>%
	ggplot(.,aes(x=log2FoldChange,y=log2(lfc_MT_MB)))+
	geom_point()+
	scale_x_continuous(paste0('mRNA Fold Change'))+
	scale_y_continuous(paste0('Silac Fold Change (Ratio.M.L.Normalized.Eluate)'))+
	ggtitle(paste0('Silac vs mRNA MT vs MB\nSignificant mRNA Change Only'),sub=cortext)+
	# geom_text(data=NULL,aes(label=cortext,x=Inf,y=-Inf),hjust=-1,vjust=1)+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))
}
}

#resultslist$MT_vs_MB                    resultslist$Ctrl_MT_A72_vs_MT0          resultslist$Aspecific_MT_A72_vs_MT0
#resultslist$Ctrl_MT_A24_vs_MT           resultslist$Aspecific_MT_A24_vs_MT      resultslist$Ctrl_MT_Amtmbsilac
#resultslist$Ctrl_MT_A72_vs_MT_A24       resultslist$Aspecific_MT_A72_vs_MT_A24
{
silaccol='Ratio.H.M.Normalized.Eluate'
mtmbsilac = msilacdata%>%select(gnm,lfc_silac = !!sym(silaccol),pval_MT_MB = Ratio.H.M.Significance.B..Eluate)
contr='FullEff_MT_A24_vs_MT'
contr='Ctrl_MT_A24_vs_MT'
mrnachagne = resultslist[[contr]]%>%select(gene_id,log2FoldChange,padj)%>%left_join(gid2gname,by=c(gene_id='g_id'))

mrna_silac_df = mtmbsilac%>%left_join(mrnachagne)

{
cortext = mrna_silac_df%>%filter(padj<0.05)%>%{cor.test(.$log2FoldChange,.$lfc_silac)}%>%tidy%>%
	{str_interp('rho = ${round(.$estimate,3)} (${round(.$conf.low,3)} - ${round(.$conf.high,3)} ),p = ${format(.$p.value,3,digits=2)}')}

#now plot
plotfile<- here(paste0('plots/','mrna_v_silac_fc_amt_mt','.pdf'))
pdf(plotfile)
print(
	mrna_silac_df%>%
	filter(padj<0.05)%>%
	ggplot(.,aes(x=log2FoldChange,y=log2(lfc_silac)))+
	geom_point()+
	scale_x_continuous(paste0('mRNA Fold Change (',contr,')'))+
	scale_y_continuous(paste0('Silac Fold Change (',silaccol,')'))+
	ggtitle(paste0('Silac vs mRNA'),sub=cortext)+
	# geom_text(data=NULL,aes(label=cortext,x=Inf,y=-Inf),hjust=-1,vjust=1)+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))
}
}
