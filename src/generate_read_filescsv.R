df = list.files(full=TRUE,rec=TRUE,patt='TS_.*fastq.gz$','../input/')%>%
	enframe('val','file')%>%
	mutate(sample=basename(dirname(dirname(dirname(file)))))

df%>%
	select(-val)%>%
	mutate(pair_id=as.numeric(as.factor(rank(str_replace(file,'_R[12]_','_')))))%>%
	mutate(pair_id=paste0(sample,'_',pair_id))%>%
	mutate(mate=str_extract(file,'(?<=_R)[12](?=_)'))%>%
	left_join(read_tsv('../ext_data/sample_names.txt')%>%
		set_colnames(c('sample','sample_id')))%>%
	select(sample_id,file,pair_id,mate)%T>%
	write_csv('read_files.csv')


rfilesampleids = read_csv('read_files.csv')%>%.$sample_id
sfilesampleids = read_csv('sample_parameter.csv')%>%.$sample_id

# sfilesampleids%>%setdiff(rfilesampleids)
# rfilesampleids%>%setdiff(sfilesampleids)