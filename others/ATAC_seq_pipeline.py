# !/usr/bin/python
# -*- coding:utf-8 -*-

import os
from multiprocessing import Pool
import time

hg38 = "/ATAC-seq/reference/hg38_UCSC/mainChr_hg38_UCSC.fa"
hg38fai = "/ATAC-seq/reference/hg38_UCSC/mainChr_hg38_UCSC.fa.fai"

workDir = "/ATAC-seq"
dataDir = "/DataStore/tgs/lung_cancer/ATAC-seq"
fastqcDir = os.path.join(workDir, "fastqc")
noadapterDir = os.path.join(workDir, "noadapter")
bamDir = os.path.join(workDir, "bam")
bedDir = os.path.join(workDir, "bed")
fragmentDir = os.path.join(workDir, "fragment")
peakDir = os.path.join(workDir, "peak")

fastqc = "/Tools/FastQC/fastqc"
cutadapt = '/venv/bin/cutadapt"
bwa = "/Tools/bin/bwa"
samtools = "/Tools/samtools1.35/samtools"
picard_markdup = "/Tools/picard-tools-1.119/MarkDuplicates.jar"
bedtools = "/Tools/bedtools2/bedtools2/bin/bedtools"
bdg2bw = "/Tools/bdg2bw"

sample_list = []


########################################################################################
# cut adapter, omit fastqc
########################################################################################
def cutadapter(storedir, sample, outdir):
    fq1 = os.path.join(storedir, sample, "%s_1.fq.gz" % sample)
    fq2 = os.path.join(storedir, sample, "%s_2.fq.gz" % sample)
    tmp_fq1 = os.path.join(outdir, "%s_noadapter_1.tmp.fq.gz" % sample)
    tmp_fq2 = os.path.join(outdir, "%s_noadapter_2.tmp.fq.gz" % sample)
    final_fq1 = os.path.join(outdir, "%s_noadapter_1.fq.gz" % sample)
    final_fq2 = os.path.join(outdir, "%s_noadapter_2.fq.gz" % sample)
    # Step 1:  cut 3' adapter from both reads
    cmd = "%s -j 2 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -q 10 --minimum-length 30 -o %s -p %s %s %s " \
          "> %s/%s_cutadapter_step1.log 2>&1" % (cutadapt, tmp_fq1, tmp_fq2, fq1, fq2, outdir, sample)
    os.system(cmd)
    # Step 2:  cut 5' adapter from both reads
    cmd = "%s -j 2 -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -q 10 --minimum-length 30 -o %s -p %s %s %s" \
          "> %s/%s_cutadapter_step2.log 2>&1" % (cutadapt, final_fq1, final_fq2, tmp_fq1, tmp_fq2, outdir, sample)
    os.system(cmd)
    os.system("rm -rf %s %s" % (tmp_fq1, tmp_fq2))


pools = Pool(2)
for sampleid in sample_list:
    sampledir = os.path.join(noadapterDir, sampleid)
    if not os.path.exists(sampledir):
        os.mkdir(sampledir)
    pools.apply_async(cutadapter, args=(dataDir, sampleid, sampledir))

pools.close()
pools.join()
del pools


########################################################################################
# bwa mapping
########################################################################################
def bwa_mapping(storedir, sample):
    ID = '%s' % sample
    PL = 'Illumina'
    CN = 'SCU'
    RG = '@RG\tID:%s\tPL:%s\tCN:%s' % (ID, PL, CN)
    fq1 = os.path.join(storedir, sample, "%s_noadapter_1.fq.gz" % sample)
    fq2 = os.path.join(storedir, sample, "%s_noadapter_2.fq.gz" % sample)
    samfile = os.path.join(bamDir, sample, "%s.sam" % sample)
    logfile = os.path.join(bamDir, sample, "%s_mapping.log" % sample)
    os.system('{bwa} mem -M -Y -t 10 -R "{RG}" {reference} {R1} {R2} 1>{out} 2>{log}'.format(
                bwa=bwa, RG=RG, reference=hg38, R1=fq1, R2=fq2, out=samfile, log=logfile))
    bamfile = os.path.join(bamDir, sample, "%s.bam" % sample)
    os.system('{samtools} view -@ 20 -S -b {sam} -o {bam}'.format(samtools=samtools, sam=samfile, bam=bamfile))
    bamfile_sort = os.path.join(bamDir, sample, "%s_sort.bam" % sample)
    os.system("{samtools} sort -@ 20 {bam} -o {bams}".format(samtools=samtools, bam=bamfile, bams=bamfile_sort))
    os.system('{samtools} index {bams}'.format(bams=bamfile_sort, samtools=samtools))
    os.system("rm %s %s" % (samfile, bamfile))


pools = Pool()
for sampleid in sample_list:
    sampledir = os.path.join(bamDir, sampleid)
    if not os.path.exists(sampledir):
        os.mkdir(sampledir)
    pools.apply_async(bwa_mapping, args=(noadapterDir, sampleid))

pools.close()
pools.join()
del pools


########################################################################################
# select bam
########################################################################################
def select_bam(sample):
    # select: deduplicate
    bam_sort = os.path.join(bamDir, sample, "%s_sort.bam" % sample)
    bam_nodup = os.path.join(bamDir, sample, "%s_nodup.bam" % sample)
    metricsfile = os.path.join(bamDir, sample, "%s_dup_metrics.txt" % sample)
    statfile = os.path.join(bamDir, sample, "%s_reads.stat" % sample)
    logfile = os.path.join(bamDir, sample, "%s_markdup.log" % sample)
    os.system("nohup java -jar {picard_markdup} INPUT={sort} OUTPUT={nodup} METRICS_FILE={metrics} \
        VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
        REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=4000000 >{log} 2>&1".format(
        picard_markdup=picard_markdup, sort=bam_sort, nodup=bam_nodup, metrics=metricsfile,log=logfile))
    os.system("{samtools} flagstat {sort} > {stat}".format(samtools=samtools, sort=bam_sort, stat=statfile))
    os.system("{samtools} flagstat {nodup} >> {stat}".format(samtools=samtools, nodup=bam_nodup, stat=statfile))
    # select: aligned;Mapping quality(MAPQ)>30;mate aligned(RNEXT);observed Template Length(TLEN)<=2000;
    # exclude multi-mapped reads(get Uniquely mapped reads) ;no chrM
    bam_select = os.path.join(bamDir, sample, "%s_select.bam" % sample)
    os.system('''%s view -F 4 -h %s|awk -F "\t" '$1~/^@/||($5>=30 && $7=="=" && $9 > -2000 && $9 < 2000){print $0}'|
    grep -v -e "XA:Z" -e "SA:Z" -e 'chrM'|%s view -Sb > %s''' % (samtools, bam_nodup, samtools, bam_select))
    os.system("{samtools} flagstat {select} >> {stat}".format(samtools=samtools, select=bam_select, stat=statfile))
    # bam to bed
    bedselect = os.path.join(bedDir, sample, "%s_select.bed" % sample)
    os.system(' {bedtools} bamtobed -i {bam} > {bed} '.format(bedtools=bedtools, bam=bam_select, bed=bedselect))


pools = Pool(2)
for sampleid in sample_list:
    sampledir = os.path.join(bedDir, sampleid)
    if not os.path.exists(sampledir):
        os.mkdir(sampledir)
    pools.apply_async(select_bam, args=(sampleid,))  # remember the ',' !!!

pools.close()
pools.join()
del pools

########################################################################################
# bed to bedpe,bedpair,fragment
########################################################################################
chromlenDic = {}
fai = open(hg38fai, "r")
linef = fai.readline()
while linef:
    contentf = linef.split()
    chromlenDic[contentf[0]] = contentf[1]
    linef = fai.readline()

fai.close()



def bed_convert(bed, bedpe, bedpair, fragment, sample):
    f = open(bed, "r")
    name_dic = {}
    line = f.readline()
    while line:
        content = line.split()
        chrom = content[0]
        start = content[1]
        end = content[2]
        name = content[3].split("/")[0]
        read_name = content[3]
        score = content[4]
        strand = content[5]
        if name not in name_dic:
            name_dic[name] = [chrom, start, end, read_name, score, strand]
        else:
            name_dic[name].extend([chrom, start, end, read_name, score, strand])
        line = f.readline()
    f.close()
    outbedpe = open(bedpe, "w")
    outbedpair = open(bedpair, "w")
    outfragment = open(fragment, "w")
    for name in name_dic:
        value = name_dic[name]
        if len(value) == 12 and value[0] == value[6]: # reads may mapping to diffenrent chrom
            newbedpe = "\t".join(value[0:3] + value[6:9] +
                                 [value[3], value[9], value[4], value[-2], value[5], value[-1]]) + "\n"
            outbedpe.write(newbedpe)
            mate1 = "\t".join(value[0:6]) + "\n"
            mate2 = "\t".join(value[6:12]) + "\n"
            outbedpair.write(mate1)
            outbedpair.write(mate2)
            info = map(int, value[1:3] + value[7:9])
            if value[5] == "+":
                info[0] += 4
                info[1] += 4
            else:
                info[0] += -5
                info[1] += -5
            if value[-1] == "+":
                info[2] += 4
                info[3] += 4
            else:
                info[2] += -5
                info[3] += -5
            temp_start = min(info)
            start = max(0, temp_start)
            temp_end = max(info)
            end = min(temp_end, int(chromlenDic[value[0]]))
            if temp_end > int(chromlenDic[value[0]]): # if? extending past the edge of the chromosome
                print("big" + "\t" + name + "\t" + sample + "\n")
            if value[0] != value[6]: # if? reads mapping to diffenrent chrom,didn't exist this time
                print("dif" + "\t" + name + "\t" + sample + "\n")
            newfrag = "\t".join([value[0], str(start), str(end), name, sample]) + "\n"
            outfragment.write(newfrag)
        else:
            newbedpe = "\t".join(value[0:3] + ['.', '-1', '-1'] +
                                 [value[3] + '.' + value[4], '.', value[5], '.']) + "\n"
            outbedpe.write(newbedpe)
    outfragment.close()
    outbedpe.close()
    outbedpair.close()
    os.system('sort -k1,1 -k2,2n %s -o %s' % (bedpe, bedpe))
    os.system('sort -k1,1 -k2,2n %s -o %s' % (fragment, fragment))
    os.system('sort -k1,1 -k2,2n %s -o %s' % (bedpair, bedpair))


pools = Pool(2)
for sampleid in sample_list:
    bed_select = os.path.join(bedDir, sampleid, "%s_select.bed" % sampleid)
    bed_pe = os.path.join(bedDir, sampleid, "%s_select.bedpe" % sampleid)
    bed_pair = os.path.join(bedDir, sampleid, "%s_select_pair.bed" % sampleid)
    sampledir = os.path.join(fragmentDir, sampleid)
    if not os.path.exists(sampledir):
        os.mkdir(sampledir)
    fragmentTxt = os.path.join(fragmentDir, sampleid, "%s_fragment.txt" % sampleid)
    pools.apply_async(bed_convert, args=(bed_select, bed_pe, bed_pair, fragmentTxt, sampleid))

pools.close()
pools.join()
del pools


########################################################################################
# bed pair call peak
########################################################################################
def callpeak_bigwig(bedpair, sample, outdir):
    os.system('nohup macs2 callpeak -f BED --nomodel --nolambda --shift 100 --extsize 200 -B --SPMR --call-summits '
              '-t {bed} --outdir {dirs} -n {sample_id} > {dirs}/{sample_id}.log 2>&1'.
              format(bed=bedpair, dirs=outdir, sample_id=sample))
    bedgraph = os.path.join(outdir, "%s_treat_pileup.bdg" % sample)
    os.system('bash %s %s %s' % (bdg2bw, bedgraph, hg38fai))


pools = Pool(2)
for sampleid in sample_list:
    bed_pair = os.path.join(bedDir, sampleid, "%s_select_pair.bed" % sampleid)
    sampledir = os.path.join(peakDir, sampleid)
    if not os.path.exists(sampledir):
        os.mkdir(sampledir)
    pools.apply_async(callpeak_bigwig, args=(bed_pair, sampleid, sampledir))

pools.close()
pools.join()
del pools
