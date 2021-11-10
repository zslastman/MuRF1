

dds = readRDS(here('data/dds.rds'))

sgroupvects = model.matrix(data=dds@colData,design(dds))%>%set_rownames(dds@colData$group)%>%{.[unique(rownames(.)),]}%>%
{m=.;lapply(1:nrow(.)%>%setNames(rownames(m)),function(i){m[i,]}) }
#
grpconvlist = list(
	 'undiff_Ctr_0'= 'MB',
	 'diff_Ctr_7'= 'MT',
	 'diff_Ctr_8'= 'MT_Ctrl24',
	 'diff_Ctr_10'= 'MT_Ctrl72',
	 'diff_atrophy_8'= 'MT_A24',
	 'diff_atrophy_10'= 'MT_A72'
)
#
dds = DESeq(dds)
#
silaccol=tmt_LogFCcols[1]
list_mrnachagne = lapply(tmt_LogFCcols,function(silaccol){
	sampgrp1 = str_extract(silaccol,'(?<=logFC\\.).*?(?=\\.over)')
	sampgrp2 = str_extract(silaccol,'(?<=\\.over\\.).*')
	list(grpconvlist[sampgrp1],grpconvlist[sampgrp2])
	#
	mrna_contrast = sgroupvects[[grpconvlist[[sampgrp2]]]] - sgroupvects[[grpconvlist[[sampgrp1]]]]
	#
	mrnares = results(dds,contrast = mrna_contrast)
	#
	mrnares%>%as.data.frame%>%rownames_to_column('g_id')%>%
		select(g_id,log2FoldChange,padj)%>%
		left_join(gid2gname)
})


filtfuns = list(all = .%>%identity%>%filter(log2FoldChange> -20,-lfc_tmt > -20) , 
	filt = .%>%filter( - lfc_tmt > -20,log2FoldChange> -20)%>%
	filter(padj<0.05 | (tmt_pval<0.05))
)
filtfun=filtfuns[[2]]
filtnm=names(filtfuns)[2]
mrna_tmt_df%>%filtfuns[[2]](.)

tmtcol = tmt_LogFCcols[[11]]

imap(filtfuns,function(filtfun,filtnm){
	lapply(tmt_LogFCcols,function(tmtcol){
		# tmtcol='logFC.undiff_Ctr_0.over.diff_Ctr_7'
		tmtpcol=tmtcol%>%str_replace('logFC','adj.P.Val')
		#
		mtmbtmt = mtmtdata%>%select(gnm,lfc_tmt = !!sym(tmtcol),tmt_pval = !!sym(tmtpcol))
		#
		mrna_tmt_df = mtmbtmt%>%tibble%>%left_join(tibble(list_mrnachagne[[tmtcol]]))
		{
		cortext = mrna_tmt_df%>%filtfun%>%
			{cor.test(.$log2FoldChange,- .$lfc_tmt)}%>%tidy%>%
			{str_interp('rho = ${round(.$estimate,3)} (${round(.$conf.low,3)} - ${round(.$conf.high,3)} ),p = ${format(.$p.value,3,digits=2)}')}
		#
		#now plot
		stopifnot(mrna_tmt_df%>%filtfun%>%.$log2FoldChange%>%min(na.rm=T)%>%`>`(-21))
		colnm = tmtcol%>%str_replace('logFC','')
		plotfile<- here(paste0('plots/','mrna_v_tmt_fc_',filtnm,'_',colnm,'.pdf'))
		pdf(plotfile)
		print(
			mrna_tmt_df%>%
			filtfun%>%
			ggplot(.,aes(x=log2FoldChange,y= - lfc_tmt))+
			geom_point(size=I(0.2))+
			scale_x_continuous(paste0('mRNA Fold Change '))+
			scale_y_continuous(paste0('- TMT Fold Change (',tmtcol,')'))+
			ggtitle(paste0('TMT vs mRNA'),sub=cortext)+
			# geom_text(data=NULL,aes(label=cortext,x=Inf,y=-Inf),hjust=-1,vjust=1)+
			theme_bw()
		)
		dev.off()
		message(normalizePath(plotfile))
		}
	})
})

