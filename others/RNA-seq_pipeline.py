import os
from multiprocessing import Pool
import time

hisat2='/Tools/hisat2-2.1.0/hisat2'
hisat2_index='/Ref/human_hg38/GRCh38/hg38_UCSC/mainChr_hg38_UCSC_tran'
stringtie='/Tools/stringtie-1.3.6.Linux_x86_64/stringtie'
gtf38='/Ref/human_hg38/GRCh38/hg38_UCSC/gencode.v33.annotation.gtf'
workdir='/TGS/RNA-seq/pipeline'
os.chdir(workdir)
sample_dict={}
ff=open('/TGS/RNA-seq/program/samplelist_20200407.txt','r')
line=ff.readline()
while line:
    l=line.strip().split('\t')
    sample_dict[l[0]]=l[1]
    line=ff.readline()
ff.close()

def pipelinerna(rawsample,newsample):
    sample_dir = workdir + "/" + newsample
    os.system('mkdir %s' % sample_dir)
    os.chdir(sample_dir)
    fastq1='/RNA-seq/X101SC19072146-Z03-F005-B4-40_20200118/CleanData/lncRNA_'+rawsample+'_1.clean.fq.gz'
    fastq2='/RNA-seq/X101SC19072146-Z03-F005-B4-40_20200118/CleanData/lncRNA_'+rawsample+'_2.clean.fq.gz'
    os.system('echo "1. hisat2  start  mapping  --  %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system('{hisat2} --dta -p 10 -q -x {hisat2_index} -1 {fastq1} -2 {fastq2} -S {sam_dir}/{sample}.sam'.format(hisat2=hisat2, hisat2_index=hisat2_index, fastq1=fastq1, fastq2=fastq2, sam_dir=sample_dir, sample=newsample))
    os.system( 'echo "1. hisat2  mapping  finish --  %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system('echo "2. sam to bam  --  %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system('samtools view -bS %s/%s.sam > %s/%s.bam' % (sample_dir, newsample, sample_dir, newsample))
    os.system('echo "3. sort bam  --  %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system('samtools sort -@ 5 %s/%s.bam -o %s/%s.sorted.bam' % (sample_dir, newsample, sample_dir, newsample))
    os.system( 'echo "4. stringtie  start  assembly  -- %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system( '{stringtie} {sample_dir}/{sample}.sorted.bam -o {sample_dir}/{sample}.out.gtf -p 10 -G {ref} -A {sample}.abund.tab '.format(stringtie=stringtie, sample_dir=sample_dir, sample=newsample, ref=gtf38))
    os.system( 'echo "4. stringtie  assembly  finish -- %s --%s" >> %s/%s.log' % (newsample, time.asctime(), sample_dir, newsample))
    os.system('rm {sam_dir}/{sample}.sam'.format(sam_dir=sample_dir,sample=newsample))
    os.system('rm %s/%s.bam'%(sample_dir,newsample))


for Rawsample in sample_dict.keys():
    Newsample=sample_dict[Rawsample]
    pipelinerna(Rawsample,Newsample)
