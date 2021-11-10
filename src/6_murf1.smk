import glob
import pandas as pd
from pathlib import Path
from functools import partial
from ipdb import set_trace


HEADIFTEST = '| head -n 400000'


def is_nonempty(file):
  assert Path(file).stat().st_size
def is_over_size(file,n):
  assert Path(file).stat().st_size > n
def newfolder(file,newroot):
  file = Path(file)
  assert snakedir in file.parents , str(snakedir) + " doesn't seem to be in parents of " + str(file)
  return str(Path(*(newroot,)+file.relative_to(snakedir).parts[1:]))

#shell.executable("/bin/bash")
# shell.prefix("set -e pipefail;")
# user set parameter


################################################################################
########load and check  pipeline configuration
################################################################################
  
configfile: "../src/config.yaml"
config['root'] = Path(config['root'])




snakedir = Path(config['root']).resolve() / 'pipeline'
assert snakedir.exists
TMPDIR = Path('../tmp')



seqfilesdf = pd.read_csv(config['sample_files'],dtype=str).set_index("sample_id", drop=False)
sampledf = pd.read_csv(config['sample_parameter']).set_index("sample_id", drop=False)

assert sampledf.sample_id.is_unique
if not 'mate' in seqfilesdf.columns: seqfilesdf['mate'] = '1'
if not 'pair_id' in seqfilesdf.columns:
    assert seqfilesdf['mate'].isin(['1']).all()
    seqfilesdf['pair_id']=seqfilesdf['sample_id']

if not 'file_id' in seqfilesdf.columns: seqfilesdf.insert(seqfilesdf.shape[1],'file_id',seqfilesdf.pair_id+'_R'+seqfilesdf.mate+'.fastq.gz',False)


assert sampledf.sample_id.is_unique
assert isinstance(seqfilesdf.file[0],str), "file column should be a string in read_files.csv"
assert not (pd.Series([Path(f).name for f in seqfilesdf.file]).duplicated().any()),"files need unique filenames"

seqfilesdf.mate = seqfilesdf.mate.fillna('1')
assert set(seqfilesdf.mate).issubset(set(['1','2']))

seqfilesdf.pair_id = seqfilesdf.pair_id.fillna(seqfilesdf.sample_id+'.fastq.gz')
seqfilesdf.file_id = seqfilesdf.file_id.fillna(seqfilesdf.sample_id+'_'+seqfilesdf.mate+'.fastq.gz')

assert seqfilesdf.file_id.is_unique,"pairid + mate combo must be unique"

print('found '+str(len(seqfilesdf.sample_id))+' sample ids')
print('for '+str(len(seqfilesdf))+' files')
n_paired = str(sum(seqfilesdf.groupby('pair_id').size() > 1))
print('of these, '+n_paired+' were paired')

for i in seqfilesdf.file: assert Path(i).stat().st_size > 100, "This file isn't of size 100 or more " +i


#For creating a pair ID column
# pairids = [Path(f).name.replace('_R2_','_R1_').replace('.fastq.gz','') for f in seqfilesdf.file]
# pairids = [re.sub('_[A-Za-z0-9]{8,20}$','',p) for p in pairids]
# #insert tinot he read file
# seqfilesdf.insert(seqfilesdf.shape[1],'pair_id', pairids, False)
#seqfilesdf.to_csv(config['sample_files'])

assert set(seqfilesdf.sample_id) == set(sampledf.sample_id), "Sample IDs need to be setequal in "+config['sample_files']+" and "+config['sample_parameter'] +": \n"+"seqfile ids " + seqfilesdf.sample_id[0:3] + "... \n" +"seqfile ids " + sampledf.sample_id[0:3] + "... "

for sample in sampledf.sample_id:
  for f in seqfilesdf.loc[[sample],'file']:

    fp = Path(f)
    assert fp.exists, f
    assert 'fastq.gz' in fp.name or 'fq.gz' in fp.name , f



samples = list(list(sampledf['sample_id'].unique()))
fastqs = list(seqfilesdf['file'].unique())
ribosamples = sampledf.sample_id[sampledf.assay=='ribo']
rnasamples = sampledf.sample_id[sampledf.assay=='total']

#for trimming CDS for riboseq
REF_orig=config['REF_orig']
GTF_orig=config['GTF_orig']


assert(Path(REF_orig).exists())
assert(Path(GTF_orig).exists()), GTF_orig + ", the GTF file, doesn't exist"


#local copies of the annotation
REF = Path(Path(re.sub(string=REF_orig,pattern='.(b)?gz$',repl=r'')).name)
GTF = Path(Path(re.sub(string=GTF_orig,pattern='.(b)?gz$',repl=r'')).name)

#REF = snakedir / Path(config['REF_orig']).stem.with_suffix('.fa').name
#GTF = snakedir / Path(config['GTF_orig']).stem.with_suffix('.gtf').name
#GFF = snakedir / Path(config['GFF_orig']).stem.with_suffix('.gff3').name
RNAFASTA = GTF.with_suffix('.fa')
CODINGFASTA=GTF.with_suffix('.coding.fa')
PROTEINFASTA=GTF.with_suffix('.protein.fa')
CDSFASTA=GTF.with_suffix('.cds.fa')
BED=GTF.with_suffix('.bed')



samplegroups = dict(sampledf['group'])

ALIGNERS = ['hisat2','star']
ALIGNER_TO_USE = config.get('aligner_to_use','star')

#
if config['TRIM_IDS']: 
  mod_id_sed_cmd = r''' sed -r 's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' '''
else:
  mod_id_sed_cmd = ' cat '



# rnasamples=[]

rule all:
  input:
    seqfilesdf.file,
    expand("{aligner}/data/{sample}/{sample}.bam", aligner = ALIGNER_TO_USE, sample = samples),
    ("multiqc/multiqc_report.html"),
    expand('salmon/data/{sample}/.done',sample=samples)
    

#things this needs - cutadapt, remove8N.pl, STAR,collapse_reads.pl,seqtk

# rule make_test_fastq:
#   input: '../gdrive/RiboSeq_Ribo_Transcriptome','../gdrive/RNASeq_Total_Transcriptome/'
#   output: 'preprocessed_reads/test/test.fastq.gz'
#   shell:r"""
#     set -e
#     mkdir -p preprocessed_reads/test/
#     seqtk sample -s100 ../gdrive/RNASeq_Total_Transcriptome/E13_1.fastq.gz 10000 | gzip > {output}
#     if [ 0 -eq $(gzip -l {output} | awk 'NR==2 {{print $2}}') ]; then rm {output} ; fi

#       """

rule copy_ref:
  input: REF_orig
  output: REF,str(REF)+'.fai'
  shell: """
      zless {REF_orig} >  {output[0]}
      samtools faidx {output[0]}
      """

rule link_in_anno:
  input: REF_orig=REF_orig,GTF=GTF,REFAI=str(REF)+'.fai'
  output: touch('annocheck.done')
  run:
    from pathlib import Path
    REFchrs = pd.read_csv(input.REFAI,sep='\t',header=None,names=['Chr','length','cumlength','V4','V5'],usecols=[0,1])
    REFchrs.in_ref = True
    #awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}'
    import subprocess
    from subprocess import Popen, PIPE
    from io import StringIO


    cmd=r"""awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}' """+str(GTF_orig)
    a=subprocess.Popen(cmd,stdout = subprocess.PIPE,shell=True)
    GTFchrs = pd.read_csv(StringIO(a.communicate()[0].decode('ascii')),sep=' ',header=None,names=['Chr','longest_in_GTF'])

    allchrs = pd.merge(REFchrs,GTFchrs,how='outer')
    allchrs['outofbounds'] = allchrs.longest_in_GTF.gt(allchrs.length) 
    allchrs['missing_in_REF'] = allchrs.length.isna()
    problems = allchrs.outofbounds | allchrs.missing_in_REF
    refonly = allchrs.longest_in_GTF.isna()
    assert not problems.any(),print("\n\n\n GTF, Chromosomes Don't quite match here: \n\n",allchrs[problems|refonly],"\n\n")

    cmd=r"""awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}' """+str(GTF_orig)
    a=subprocess.Popen(cmd,stdout = subprocess.PIPE,shell=True)
    GTFchrs = pd.read_csv(StringIO(a.communicate()[0].decode('ascii')),sep=' ',header=None,names=['Chr','longest_in_GTF'])

    allchrs = pd.merge(REFchrs,GTFchrs,how='outer')
    allchrs['outofbounds'] = allchrs.longest_in_GTF.gt(allchrs.length) 
    allchrs['missing_in_REF'] = allchrs.length.isna()
    problems = allchrs.outofbounds | allchrs.missing_in_REF
    assert not problems.any(),print("\n\n\n GTF, Chromosomes Don't quite match here: \n\n",allchrs[problems],"\n\n")


def name_preprocessed_reads(wc):
  filedf = seqfilesdf[seqfilesdf.sample_id==wc['sample']]  
  assert (filedf.file_id  == wc['fastq']).sum()==1, "Fastq file isn't amongst the file ids"
  filedf = filedf[filedf.file_id==wc['fastq']]

  return filedf.file

rule link_in_files:
  input: name_preprocessed_reads
  output: 'preprocessed_reads/{sample}/{fastq}'
  run:  
    sample = wildcards['sample']
    fastq = wildcards['fastq']
    shell(r"""
      mkdir -p $(dirname {output})
      ln -sf $(readlink -f {input} {output} )
    """)

#this rule is the 'signal spliter where we go from sample to indiv fastqs
def choose_processed_reads(wc,config=config):
  #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  filedf = filedf[filedf.file_id==wc['fileid']]
  isrna = not 'ribo' in sampledf.loc[wc['sample'],'assay']
  if isrna:
    out =  [config['FILT_RNA_FOLDER']+'/'+wc['sample']+'/'+f for f in filedf.file_id]
  else:
    out =  [config['FILT_RIBO_FOLDER']+'/'+wc['sample']+'/'+f for f in filedf.file_id]
  # print(wc)
  # print(out)

  return out

rule link_processed_reads:
  # input: choose_processed_reads
  # input: choose_processed_reads
  input: name_preprocessed_reads
  output: 'processed_reads/{sample}/{fastq,\w.*}'
  run: 
    shell(r"""
        mkdir -p $(dirname {output} )
        for i in $(readlink -f {input} ); do  ln -rifs $i  {output} ; done
    """)
    assert Path(output[0]).stat().st_size > 100

################################################################################
########Annotation
################################################################################
  

rule makeGTF:
  input: GTF=GTF_orig
  output: GTF
  conda: '../envs/gffread'
  #conda: '~/miniconda3/envs/seq/bin/gffread'
  shell: r""" 
      # set -x
      #with filtering output all sequences
      zless {input.GTF} \
      | {mod_id_sed_cmd} \
      | gffread -F -T -o {GTF}

    """

# rule makeCDSGTF:
#   input: REF=REF_orig,GTF=GTF_orig
#   output: CDSGTF,RNAFASTA,CDSFASTA,BED
#   conda: '../envs/gffread'
#   #conda: '~/miniconda3/envs/seq/bin/gffread'
#   shell: r""" 
#       #needs gff - output exon sequences
#       #conda activate seq
#      zless {input.GTF} | grep -P -e'\texon\t|^\#' | gffread - -F -E -g {REF} -W -w {RNAFASTA} -o /dev/null
#       #Note we are now minus the transcript and exon entries for these
#       #now make GTF
# #
#       #| grep -P -e'\tCDS\t|^\#' 
#      #with filtering, output the coding sequences filteirng out the ones that aren't in frame, have a stop codon, are pseudogenes etc.
# #      
#       zless {input.GTF} \
#         | {mod_id_sed_cmd} \
#         | gffread - -C -V -J --no-pseudo -F -E -g {REF} \
#         -W -w {CODINGFASTA} -x {CDSFASTA} -y {PROTEINFASTA} -T \
#         -o /dev/stdout > {CDSGTF}
# #
#       #now make bed
#       zless {input.GTF} | awk '{{print $1,$4,$5,"name",$6,$7}}' > {BED}
#     """


 
rule make_utrs:
  input: GTF=GTF_orig
  output: fputrs='fputrs.gtf',tputrs='tputrs.gtf'
  # script: 'make_utrfiles.R'
  run:
    shell(r"""
      set -ex
      #with filtering output all sequences
      cat {input.GTF}  \
      | awk -v OFS="\t"  '{{if($3=="five_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.fputrs} 

      cat {input.GTF} \
      | awk -v OFS="\t"  '{{if($3=="three_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.tputrs}

     
      """) 


################################################################################
########Quality checking
################################################################################
  
def get_sample_fastqs(wc,mate='1',folder='processed_reads',seqfilesdf=seqfilesdf):
   #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  filedf = filedf.loc[filedf.mate==mate,]
  # print('getting fastqs')
  # print(filedf)
  matefiles = [folder+'/'+wc['sample']+'/'+f for f in filedf.file_id]

  return(matefiles)

get_sample_fastqs2 = partial(get_sample_fastqs,mate='2')

#TODO
rule fastqc:
     input:
        lfastqs=get_sample_fastqs,
        rfastqs=get_sample_fastqs2,
     output: touch('fastqc/data/{sample}/.done')
     threads: 4
     log:'fastqc/reports/{sample}/fastqc.log'
     params:
        outdir = lambda wc: 'fastqc/data/'+wc.sample+'/'
     shell: '''
          OUTDIR=$(dirname {output[0]})
          mkdir -p {params.outdir}
          wait $(for i in {input.lfastqs} {input.rfastqs}; do $( fastqc -o {params.outdir} $i ) & done) 
        '''

rule collect_fastqc:
     input:
          all_results = expand("fastqc/data/{sample}/.done", sample=samples)
     output:
          result='fastqc/summary/fastqc_summary.tsv',
          log='fastqc/summary/fastqc_summary.log'
     shell:
          r"""
          set -e
          mkdir -p $(dirname {output.result}) 
          {SCRIPTDIR}/collect_fastqc_results.sh -i fastqc/ \
          > {output.result} \
          2> {output.log} 
          """

#note that it soft clips by default
#k -number of multimaps reported  


################################################################################
########STAR
################################################################################
  
rule star_index:
  input: REF=REF,GTF=GTF
  output: touch('starindex/.done')
  # conda:'../envs/star'
  threads: 8
  shell:r"""
      STAR \
      --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir $(dirname {output}) \
      --sjdbGTFfile {input.GTF} \
      --genomeFastaFiles {input.REF}
      """


rule star:
     input:
          lfastqs=get_sample_fastqs,
          rfastqs=get_sample_fastqs2,
          STARINDEX='starindex/.done',
          # bowtie_index='bowtie_index/.done',
     output:
          done = touch('star/data/{sample,[^/]+}/.done'),bam='star/data/{sample}/{sample}.bam'
     threads: 8
     run:
          input.STARINDEX=input.STARINDEX.replace('.done','')
          markdup = '' if sampledf.assay[wildcards['sample']] == 'ribo' else '-m'
          platform = 'NotSpecified'
          outputdir = os.path.dirname(output[0])
          read_pattern = sampledf.read_pattern[wildcards['sample']]
          
          # def get_file_string(wc,seqfilesdf=seqfilesdf)
          # #correct zcat strings based on read pairs
          filedf = (seqfilesdf.loc[[wildcards['sample']]])

          lfilestring = '<(zcat '+' '.join(input.lfastqs)+')' 
          rfilestring = '<(zcat '+' '.join(input.rfastqs)+')'  if input.rfastqs   else '' 

          filestring = lfilestring+' '+rfilestring


          # print(lfilestring)
          # print(rfilestring)
          # print(filestring)

          assert(len(input.lfastqs)>1)
          # assert(len(input.lfastqs)>1)

          repdir = outputdir.replace('data','reports')
          # tophatindex =input['bowtie_index'].replace('.done','')
          
          halfthreads = threads/2
          sortmem = str(int(5000/halfthreads))+'M'

          # remap = '1' if sampledf.assay[wildcards['sample']] == 'ribo' else ''
          remap = '' 
          
          sample = wildcards['sample']
          shell(r"""
            set -xeo pipefail
         MY_TMP_DIR=$(mktemp -d)

        # trap "set -x; rm -rf ${{MY_TMP_DIR}}" EXIT KILL TERM INT HUP

         mkdir -p $MY_TMP_DIR
        mkdir -p $MY_TMP_DIR/star
        mkdir -p $MY_TMP_DIR/tophat2

        #--outSAMmultNmax 20 --winAnchorMultimapNmax 50 --outFilterMultimapNmax 20 \
        #--genomeLoad NoSharedMemory \
        #  --limitOutSJcollapsed 20000000 \
        #  --limitOutSJoneRead 20000000 \
               
        STAR \
              --genomeDir $(realpath {input.STARINDEX})  \
              --runThreadN {threads} \
              --outSAMunmapped Within \
              --outFilterType BySJout \
              --outMultimapperOrder Random \
              --alignSJoverhangMin 8 \
              --alignSJDBoverhangMin 1 \
              --outFilterMismatchNmax 999 \
              --outFilterMismatchNoverLmax 0.04 \
              --alignIntronMin 20 \
              --alignIntronMax 5000000 \
              --alignMatesGapMax 5000000 \
              --quantMode GeneCounts \
              --outSAMattributes NH HI AS NM MD \
              --outSAMtype BAM  Unsorted\
              --outSAMattrRGline \"ID:{sample}\" \"SM:{sample}\" \"PL:{platform}\" \
              --outFileNamePrefix ${{MY_TMP_DIR}}/star/ \
              --outReadsUnmapped Fastx \
              --readFilesIn {filestring}
          
                
         samtools sort \
          -@ {halfthreads}\
          -m {sortmem} \
          -T ${{MY_TMP_DIR}} \
          -o {outputdir}/{sample}.bam \
          ${{MY_TMP_DIR}}/star/Aligned.out.bam
      
        samtools index {outputdir}/{sample}.bam 

        mkdir -p {repdir}
        samtools stats {outputdir}/{sample}.bam > {repdir}/{sample}.bamstats.txt
        samtools flagstat {outputdir}/{sample}.bam > {repdir}/{sample}.flagstat.log
        samtools idxstats {outputdir}/{sample}.bam > {repdir}/{sample}.idxstats.log
        
        cp  ${{MY_TMP_DIR}}/star/ReadsPerGene.out.tab {outputdir}/ReadsPerGene.out.tab
        cp  ${{MY_TMP_DIR}}/star/SJ.out.tab {outputdir}/
        cp  ${{MY_TMP_DIR}}/star/{{Log.final.out,Log.out}} {repdir}/

          """)


def get_multiqc_dirs(wildcards,input):
      reportsdirs = list(input)
      reportsdirs=[s.replace(ALIGNER_TO_USE+'/data',ALIGNER_TO_USE+'/reports') for s in reportsdirs]
      reportsdirs=[s.replace('tophat2/data','tophat2/reports') for s in reportsdirs]
      reportsdirs=[os.path.dirname(s) for s in list(reportsdirs)]
      assert len(reportsdirs) > 0
      return(reportsdirs)

rule multiqc:
  input:
      expand("fastqc/data/{sample}/.done", sample = samples),
      expand(ALIGNER_TO_USE+"/data/{sample}/{sample}.bam", sample = samples),
      expand("qc/data/{sample}/.done", sample = samples),
      # expand("tophat2/data/{sample}/.done", sample = samples),
      # [f.replace('input','filter_reads') for f in  seqfilesdf.file[ribosamples]],
      # expand("feature_counts/data/{sample}/feature_counts", sample = samples),


      # 'sample_file.txt'
  conda: '../envs/multiqc'
  params: 
    multiqcscript = config['multiqcscript'],
    sample_reads_file=config['sample_files'],
    reportsdirs= get_multiqc_dirs,
    sampnames = '--sample-names '+config.get('samplenamesfile') if config.get('samplenamesfile') else ''
  output:
    'multiqc/multiqc_report.html'
  shell:r"""
      cat {params.sample_reads_file} | sed 's/.fastq.gz//g' | sed 's/\t.*\//\t/g' \
      | awk -vOFS='\t' 'BEGIN{{print "fastqname","samplename"}}{{sumsamp[$1] = sumsamp[$1]+1;print $2,$1"_fq"sumsamp[$1]}}' \
      > multiqc/samplenames.txt

      {params.multiqcscript} {params.reportsdirs} -fo $(dirname {output[0]}) {params.sampnames}
      """

################################################################################
########ISoform specific quant
################################################################################
  
TRFASTA='../ext_data/gencode.vM12.transcripts.fa'

rule make_salmon_index:
  threads: 8
  input: TRFASTA
  output: salmonindex = touch('salmonindex/.done')
  run:
    shell(""" salmon index -p {threads} --perfectHash -k 21 -t {input} -i salmonindex""")

SALMONLIBTYPE = config['SALMONLIBTYPE']
rule salmon:
  input:
    lfastqs=get_sample_fastqs,
    rfastqs=get_sample_fastqs2,
    salmonindex='salmonindex/.done'
  params:
    salmonindex = lambda wc,input: input.salmonindex.replace('.done',''),
  output:
      done = touch('salmon/data/{sample}/.done')
  threads: 4
  shell:r"""
      set -ex
      mkdir -p salmon/reports/{wildcards.sample}
      mkdir -p salmon/data/{wildcards.sample}
      salmon quant \
      -p {threads} \
      -l {SALMONLIBTYPE} \
      --seqBias \
      --gcBias \
      -i {params.salmonindex} \
      -1 <(zcat {input.lfastqs} ) -2 <(zcat {input.rfastqs})\
      --output salmon/data/{wildcards.sample} \
      --validateMappings
"""

rule star_transcript:
  input:
    fastq='processed_reads/{sample}/{sample}.fastq.gz',
    transcriptindexfold='deepshape/StarIndex/transcript',
  output: 
    bam='star_transcript/data/{sample}/{sample}.sort.bam',
    bai='star_transcript/data/{sample}/{sample}.sort.bam.bai',
  params:
    bamnosort=lambda wc,output: output.bam.replace('sort.','')
  shell:r"""

  mkdir -p star_transcript/data/{wildcards.sample}/{wildcards.sample}.fastq.gz

  STAR --runThreadN 15 --genomeDir {input.transcriptindexfold} \
    --readFilesIn <(zcat {input.fastq}) \
    --outFileNamePrefix star_transcript/{wildcards.sample}.transcript_ \
    --outSAMtype BAM Unsorted \
    --outSAMmode NoQS \
    --outSAMattributes NH NM \
    --seedSearchLmax 10 \
    --outFilterMultimapScoreRange 0 \
    --outFilterMultimapNmax 255 \
    --outFilterMismatchNmax 1 \
    --outFilterIntronMotifs RemoveNoncanonical

  mv star_transcript/{wildcards.sample}.transcript_Aligned.out.bam {params.bamnosort}
  samtools sort {params.bamnosort} -o {output.bam}
  samtools index {output.bam}

  """
