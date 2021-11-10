# MuRF1 Project Git

This is the codebase for my collaboration with Suzuka Nagawa on her study of MuRF1 function in muscular atrophy. The project consists of a standard RNAseq expression analysis and integration with proteomic data. The files are listed below, in the order in which they are run.:

used by various scripts
``` src/Rprofile.R```

this is the snakemake pipeline for the basic analysis, does alignment,
transcript quantification (salmon) and some QC
```
src/0_6_murf1.smk  
uses:
	 src/sample_parameter.csv
	 src/snake_job.sh
 	 src/config.yaml
	 src/config_pipeline.json
 ```




This runs deseq, carries out various things like QC and go term analyses
creates objects used by other scripts (is isoform aware)
```
src/4_rnaseq_de_report_2.Rmd
uses:
	 src/data/GTOGO.rds
 	src/genesofinterest.txt
 	src/gofuncs.R
```
this carries out rnaseq analysis but isoform *specific* (not just aware)
```src/4b_rnaseq_dtu_report.Rmd```


plot expression of specific scripts.
 ```src/3b_exprplots.R```

create plots of ms and rnaseq log fold changes
 ```src/3_plot_ms_mrna_scatters.r```

this uses proteomics and transcirptoimcs to produce pltos of specific loci
```src/5b_compressed_locusplot.R```

load the ms data in final normalized form for use with quadrant analysis
``` src/2_proteome_normalization.R```

carry out linear model of ms and rnaseq, define quadrants
```src/1_quadrant_analysis.R```

this plots the quadrant data in an interactive way.
``` src/0_plotly.Rmd```

this exports various files to the project dropbox (mounted locally) 
 ```src/populate_dropbox.sh```


 Notes for likely future modification.
 The 'proteome_normalisation' script is likely to be the one you want to modify if you want to change e.g. proteomics normalisation.
 The 'quandrant_analysis script' (search for 'case_when') is where changes to the quadrant definitions will go.