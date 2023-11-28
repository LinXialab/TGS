import os
import pandas as pd
from multiprocessing import Pool
import time
import pysam
import re
import numpy as np
import cairosvg
import genomeview
from collections import Counter

hg38_resultDir = '/rna/RNA-pipeline_bam'
ins_locdir = '/Assembly_result'
dup_locdir = '/Assembly_result'

bedtools = '/biosoft/bedtools-2.28.0/bin/bedtools'
samtools = '/biosoft/samtools-1.9/samtools'
hisat2 = '/software/hisat2-2.1.0/hisat2'
hisat2_build = '/software/hisat2-2.1.0/hisat2-build'
stringtie = '/software/stringtie-2.1.4/stringtie'
######################################################################################################

# the site of TE on INS-assembled sequence (bed)
tf_te_info = pd.read_csv("/INSandDUP_V3/INS_11159/11159_INS_reannotate_duplication_repeatinfo.tsv",
                         sep='\t',header=None)
tf_te_info['type'] = ['DUP']*tf_te_info.shape[0]
tf_te_info1 = pd.read_csv("/INSandDUP_V3/INS_11159/11159_INS_reannotate_insertion_repeatinfo.tsv",
                         sep='\t',header=None)
tf_te_info1['type'] = ['INS']*tf_te_info1.shape[0]
tf_te_info = tf_te_info.append(tf_te_info1)
out = open("all.INS.DUP.anno.region.clean.te_list.bed",'w')
m=0
for i in range(0, tf_te_info.shape[0]):
    info = tf_te_info.iloc[i].map(str).values.tolist()
    refname = info[0]
    sample = info[0].split("_")[0]
    type = info[-1]
    te_detail = info[3]
    if '|' in teinfo1[teinfo1['newid_xl']==refname].iloc[0].values.tolist()[2]:
        te_type = 'complex'
        m+=1
    else:
        te_type = info[1]
    if te_detail == 'N':
        te_detail = info[1]
    locfile = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
    if not os.path.exists(locfile):
        locfile = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
    new_ins_region = open(locfile).readline().strip().split("\t")
    new_ins_region[1] = str(int(new_ins_region[1])-50)
    new_ins_region[2] = str(int(new_ins_region[2])+50)
    te = info[2]
    if te != 'Unmask':
        if int(new_ins_region[1])+int(te.split('_')[1])-1<=int(new_ins_region[2]):
            newline = '\t'.join([refname, str(int(new_ins_region[1])+int(te.split('_')[0])-1),
                                str(int(new_ins_region[1])+int(te.split('_')[1])-1),type,sample,te_type,
                                info[1],te_detail])+'\n'
        else:
            newline = '\t'.join([refname, str(int(new_ins_region[1]) + int(te.split('_')[0]) - 1),
                                str(int(new_ins_region[1]) + int(te.split('_')[1]) - 1), type, sample,te_type,
                                info[1],te_detail]) + '\n'
    else:
        newline = '\t'.join([refname, new_ins_region[1],
                            new_ins_region[2], type, sample,te_type,
                            info[1],te_detail]) + '\n'
    out.write(newline)


out.close()

######################################################################################################


######################################################################################################
sample_list = {}
for Sample in os.listdir("/rna/extract_hg38/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/rna/extract_hg38/' + Sample
            sample_list[Sample] = sampledir

for Sample in os.listdir("/rna/extract_hg38_1/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/rna/extract_hg38_1/' + Sample
            sample_list[Sample] = sampledir

for Sample in os.listdir("/rna/extract_hg38_2/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/rna/extract_hg38_2/' + Sample
            sample_list[Sample] = sampledir




######### all RNA reads
hg38_te = 'hg38.repeatmasker.reannot.sorted.bed'
rna_fq_dir = '/nanopore_wxy/RNA-seq/cleandata'
workdir = '/nanopore/xl/INS_DUP_RNA'
os.chdir(workdir)


## -- extracting split reads from all reads
def all_split_reads(sample, bam):
    samfile = pysam.AlignmentFile(bam, 'rb')
    outbam = pysam.AlignmentFile(sample + '_allreads_split.bam', 'wb', template=samfile)
    allreads = samfile.fetch()
    split_reads = open(sample + "_split_readsname.txt", "w")
    for read in allreads:
        if not read.is_unmapped:
            cigar = read.cigarstring
            read_name = read.query_name
            if read.is_read1:
                read_name = read_name + '/1'
            else:
                read_name = read_name + '/2'
            if 'N' in cigar:
                split_reads.write(read_name + '\t%s\n' % cigar)
                outbam.write(read)
    samfile.close()
    outbam.close()
    os.system('sort -u %s_split_readsname.txt > %s_split_readsname_uniq.txt' % (sample, sample))
    os.system(bedtools + ' bamtobed -split -i %s_allreads_split.bam > %s_allreads_split.bed' % (sample, sample))


## -- mapping the discondant reads on INS sequence
def split_reads(sample, bam):
    log = open(sample + '_split.log', 'w')
    log.write(sample + " " + time.asctime() + ' start\n')
    os.system(samtools + ' index ' + bam)
    samfile = pysam.AlignmentFile(bam, 'rb')
    outbam = pysam.AlignmentFile(sample + '_split_disc.bam', 'wb', template=samfile)
    allreads = samfile.fetch()
    split_reads = open(sample + "_split.disc_readsname.txt", "w")
    for read in allreads:
        if not read.is_unmapped:
            cigar = read.cigarstring
            read_name = read.query_name
            if read.is_read1:
                read_name = read_name + '/1'
            else:
                read_name = read_name + '/2'
            if 'N' in cigar:
                split_reads.write(read_name + '\t%s\n' % cigar)
                outbam.write(read)
    samfile.close()
    outbam.close()
    split_reads.close()
    log.write(sample + " " + time.asctime() + ' end\n')
    log.close()
    os.system('sort -u %s_split.disc_readsname.txt > sorted_uniq_%s_split.disc_readsname.txt' % (
        sample, sample))


# 确定split-reads拆分的segment比对到哪些位置
def segment(sample):
    split_read1 = pd.read_csv(os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_split.disc_readsname.txt'),
                              sep='\t', header=None, names=['read', 'cigar'])
    reads_ins = pd.read_csv(os.path.join(workdir, sample, sample + '_ins_disc.reads.txt'),
                            sep='\t', header=None, names=['ins', 'start', 'end', 'read', 'MAPQ', 'strand'])
    new_reads_ins = reads_ins.copy()
    delREAD = []
    split_bed = open(os.path.join(workdir, sample, sample + '_ins_split-reads.disc_detail.bed'), 'w')
    reads = list(split_read1['read'])
    for read in reads:
        # for indx in ['/1', '/2']:
        sub_split_reads_ins = reads_ins[reads_ins['read'] == read]
        if len(sub_split_reads_ins) != 0:
            delREAD.append(read)
            chr = sub_split_reads_ins.iloc[0, 0]
            start = sub_split_reads_ins.iloc[0, 1]
            readid = read
            sub_split = split_read1[split_read1['read'] == readid]
            cigar_list = re.split(r'(\d+)', str(list(sub_split.cigar)[0]))
            del cigar_list[0]
            cigar_list1 = [cigar_list[i:i + 2] for i in range(0, len(cigar_list), 2)]
            for n in range(0, len(cigar_list1)):
                if 'M' in cigar_list1[n]:
                    end = start + int(cigar_list1[n][0])
                    newline = '\t'.join([chr, str(start), str(end), readid]) + '\n'
                    split_bed.write(newline)
                elif 'N' in cigar_list1[n]:
                    start = end + int(cigar_list1[n][0])
    split_bed.close()
    delREAD = np.unique(delREAD)
    new_reads_ins = new_reads_ins[~new_reads_ins['read'].isin(delREAD)]
    new_reads_ins.to_csv(os.path.join(workdir, sample, sample + '_ins_reads.disc_nonsplit.txt'), sep='\t', header=None,
                         index=None)
    os.system('sort -u %s | sort -k1,1 -k2,2n > %s' % (
    os.path.join(workdir, sample, sample + '_ins_split-reads.disc_detail.bed'),
    os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_ins_split-reads.disc_detail.bed')))


def new_mapping(sample):
    fq1 = sample + '.disc.Read1.fastq'
    fq2 = sample + '.disc.Read2.fastq'
    new_subworkdir = ' /new_somatic_SV_TGS_bed/INS_DUP_RNA/%s' % sample
    sam = sample + '_discReads_mapped.sam'
    os.system(hisat2 + ' -p 10 -x {sample_workdir} -1 {fq1} -2 {fq2} -S {sam}'.format(
        sample_workdir=new_subworkdir, fq1=fq1, fq2=fq2, sam=sam))
    os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (sam, sam.split(".sam")[0]))
    os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (sam.split(".sam")[0], sam.split(".sam")[0]))
    # os.system('rm *sam')
    split_reads(sample, 'sorted_%s.bam' % sam.split(".sam")[0])  # 提取bam中split reads
    os.system(bedtools + ' bamtobed -i sorted_{idx}.bam > sorted_{idx}.bed'.format(idx=sam.split(".sam")[0]))
    os.system(
        bedtools + ' intersect -a sorted_{idx}.bed -b {dir}/{sample}_ins_v2n.bed -wa | sort -u > {sample}_ins_disc.reads.txt'.format(
            idx=sam.split(".sam")[0], sample=sample, dir=new_subworkdir))
    segment(sample)  # ins中split reads比对到组装INS什么位置
    os.system(bedtools + ' intersect -a sorted_uniq_{sample}_ins_split-reads.disc_detail.bed \
                    -b {dir}/{sample}_ins_v2n.bed -wa | sort -u > \
                    {sample}_ins_split-reads.disc_detail.txt'.format(sample=sample,
                                                                     dir=new_subworkdir))  # 有segment落在INS内的split-read


def preparation(sample, bam):
    subworkdir = os.path.join(workdir, sample)
    if not os.path.exists(subworkdir):
        os.mkdir(subworkdir)
    os.chdir(subworkdir)
    # all_split_reads(sample, bam)
    os.system(samtools + ' sort -n -@ 10 %s > %s_RNA_sorted_name.bam' % (bam, sample))
    os.system(samtools + ' view -h {sample}_RNA_sorted_name.bam| samblaster -e -o {sample}_tmp.sam \
                            -d {sample}.disc.sam -s {sample}.split.sam'.format(
        sample=sample))
    os.system(samtools + ' view -@ 10 -bS {sample}.disc.sam | {samtools} sort -n -@ 10 > {sample}.disc.bam'.format(
        sample=sample, samtools=samtools))
    os.system('rm *.sam')
    os.system(
        bedtools + ' bamtofastq -i {sample}.disc.bam -fq {sample}.disc.Read1.fastq -fq2 {sample}.disc.Read2.fastq'.format(
            sample=sample))
    split_bed = sample + '_allreads_split.bed'
    os.system(bedtools + ' intersect -a {hg38_te} -b {split_bed} -wa -wb > {sample}_split-reads_hg38te.txt'.format(
        sample=sample, split_bed=split_bed, hg38_te=hg38_te))
    new_mapping(sample)


pools = Pool(5)
for sample in list(sample_list.keys()):
    bam = '/rna/RNA-pipeline_bam/%s.sorted.bam' % sample
    pools.apply_async(preparation, args=(sample, bam,))

pools.close()
pools.join()
del pools

seqtk = '/tools/seqtk/seqtk'
workdir = '/INS_DUP_RNA'
exonBed = '/ref/sorted_gencode.v33.annotation_exon.bed'
rna_fq_dir = '/nanopore_wxy/RNA-seq/cleandata'


# extract split reads
## -- mapping the discondant reads on INS sequence
def new_split_reads(sample, bam):
    log = open(sample + '_new_splitreads.log', 'w')
    log.write(sample + " " + time.asctime() + ' start\n')
    os.system(samtools + ' index ' + bam)
    samfile = pysam.AlignmentFile(bam, 'rb')
    outbam = pysam.AlignmentFile(sample + '_split_all-split.bam', 'wb', template=samfile)
    allreads = samfile.fetch()
    split_reads = open(sample + "_split.all-split_readsname.txt", "w")
    for read in allreads:
        if not read.is_unmapped:
            cigar = read.cigarstring
            read_name = read.query_name
            if read.is_read1:
                read_name = read_name + '/1'
            else:
                read_name = read_name + '/2'
            if 'N' in cigar:
                split_reads.write(read_name + '\t%s\n' % cigar)
                outbam.write(read)
    samfile.close()
    outbam.close()
    split_reads.close()
    log.write(sample + " " + time.asctime() + ' end\n')
    log.close()
    os.system('sort -u %s_split.all-split_readsname.txt > sorted_uniq_%s_split.all-split_readsname.txt' % (
        sample, sample))


# 确定split-reads拆分的segment比对到哪些位置
def new_segment(sample):
    split_read1 = pd.read_csv(os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_split.all-split_readsname.txt'),
                              sep='\t', header=None, names=['read', 'cigar'])
    reads_ins = pd.read_csv(os.path.join(workdir, sample, sample + '_ins_all-split.reads.txt'),
                            sep='\t', header=None, names=['ins', 'start', 'end', 'read', 'MAPQ', 'strand'])
    new_reads_ins = reads_ins.copy()
    delREAD = []
    split_bed = open(os.path.join(workdir, sample, sample + '_ins_split-reads.all-split_detail.bed'), 'w')
    reads = list(split_read1['read'])
    for read in reads:
        # for indx in ['/1', '/2']:
        sub_split_reads_ins = reads_ins[reads_ins['read'] == read]
        if len(sub_split_reads_ins) != 0:
            delREAD.append(read)
            chr = sub_split_reads_ins.iloc[0, 0]
            start = sub_split_reads_ins.iloc[0, 1]
            readid = read
            sub_split = split_read1[split_read1['read'] == readid]
            cigar_list = re.split(r'(\d+)', str(list(sub_split.cigar)[0]))
            del cigar_list[0]
            cigar_list1 = [cigar_list[i:i + 2] for i in range(0, len(cigar_list), 2)]
            for n in range(0, len(cigar_list1)):
                if 'M' in cigar_list1[n]:
                    end = start + int(cigar_list1[n][0])
                    newline = '\t'.join([chr, str(start), str(end), readid]) + '\n'
                    split_bed.write(newline)
                elif 'N' in cigar_list1[n]:
                    start = end + int(cigar_list1[n][0])
    split_bed.close()
    delREAD = np.unique(delREAD)
    new_reads_ins = new_reads_ins[~new_reads_ins['read'].isin(delREAD)]
    new_reads_ins.to_csv(os.path.join(workdir, sample, sample + '_ins_reads.all-split_nonsplit.txt'), sep='\t',
                         header=None, index=None)
    os.system('sort -u %s | sort -k1,1 -k2,2n > %s' % (
    os.path.join(workdir, sample, sample + '_ins_split-reads.all-split_detail.bed'),
    os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_ins_split-reads.all-split_detail.bed')))


def new_mapping(sample, fq1, fq2):
    new_subworkdir = ' /new_somatic_SV_TGS_bed/INS_DUP_RNA/%s' % sample
    sam = sample + '_splitReads_mapped.sam'
    os.system(hisat2 + ' -p 10 -x {sample_workdir} -1 {fq1} -2 {fq2} -S {sam}'.format(
        sample_workdir=new_subworkdir, fq1=fq1, fq2=fq2, sam=sam))
    os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (sam, sam.split(".sam")[0]))
    os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (sam.split(".sam")[0], sam.split(".sam")[0]))
    # os.system('rm *sam')
    new_split_reads(sample, 'sorted_%s.bam' % sam.split(".sam")[0])  # 提取bam中split reads
    os.system(bedtools + ' bamtobed -i sorted_{idx}.bam > sorted_{idx}.bed'.format(idx=sam.split(".sam")[0]))
    os.system(
        bedtools + ' intersect -a sorted_{idx}.bed -b {dir}/{sample}_ins_v2n.bed -wa | sort -u \
                    > {sample}_ins_all-split.reads.txt'.format(
            idx=sam.split(".sam")[0], sample=sample, dir=new_subworkdir))
    new_segment(sample)  # ins中split reads比对到组装INS什么位置
    os.system(bedtools + ' intersect -a sorted_uniq_{sample}_ins_split-reads.all-split_detail.bed \
                    -b {dir}/{sample}_ins_v2n.bed -wa | sort -u > \
                    {sample}_ins_split-reads.all-split_detail.txt'.format(sample=sample,
                                                                          dir=new_subworkdir))  # 有segment落在INS内的split-read


def extract_mapped_split_reads(sample):
    subworkdir = os.path.join(workdir, sample)
    if not os.path.exists(subworkdir):
        os.mkdir(subworkdir)
    os.chdir(subworkdir)
    splitbed = sample + '_allreads_split.bed'
    nonExon_splitbed = 'nonExon' + sample + '_allreads_split.bed'
    os.system('sort -k1,1 -k2,2n %s|' % splitbed + bedtools + ' intersect -a - -b %s -v -sorted > %s' % (
    exonBed, nonExon_splitbed))
    namelst = '%s_split_readsname_uniq.lst' % sample
    os.system("awk -F '/' '{print $1}' %s | awk -F '\t' '{print $4}' | sort -u > %s" % (nonExon_splitbed, namelst))
    fq1 = os.path.join(rna_fq_dir, sample, sample, sample + '_1.fq.gz')
    fq2 = os.path.join(rna_fq_dir, sample, sample, sample + '_2.fq.gz')
    fq_list = {'fq1': fq1, 'fq2': fq2}
    for fq in fq_list:
        outfq = fq_list[fq].split("/")[-1].split('.fq')[0] + '_splitReads.fq'
        os.system('gunzip -c %s | ' % fq_list[fq] + seqtk + ' subseq - %s > %s' % (namelst, outfq))
    new_mapping(sample, sample + '_1_splitReads.fq', sample + '_2_splitReads.fq')


pools = Pool(5)
for sample in list(sample_list.keys()):
    pools.apply_async(extract_mapped_split_reads, args=(sample,))

pools.close()
pools.join()
del pools

# 统计INS内RNA-seq reads数量
genebed = '/ref/sorted_gencode.v33.annotation_gene_20200312.bed'
gene5kb_bed = '/ref/gencode_v33_annotation_gtf_gene5kb_cosmic-ncg6.bed'

for idx in ['all-split', 'disc']:
    reads_result = open('samples_reads_INS_stat_%s.txt' % idx, 'w')
    reads_result.write('sample\tchr\tstart\tend\tid\treads\tsplit_reads\n')
    for sample in list(sample_list):
        if ('Lung' in sample) and ("ht2" not in sample):
            # os.system(samtools + ' flagstat {workdir}/{sample}/sorted_{sample}_unmapReads_mapped.bam > \
            # {workdir}/{sample}/sorted_{sample}_unmapReads_mapped.flagstat'.format(
            #    sample=sample,workdir=workdir))
            reads = os.path.join(workdir, sample,
                                 sample + '_ins_reads.%s_nonsplit.txt' % idx)
            v2nBed_break = os.path.join(workdir, sample,
                                        sample + '_ins_v2n_break.bed')
            newReads = os.path.join(workdir, sample,
                                    sample + '_ins_reads.%s_nonsplit_in_INS.txt' % idx)
            split = os.path.join(workdir, sample,
                                 sample + '_ins_split-reads.%s_detail.txt' % idx)
            if os.path.exists(reads) and os.path.getsize(reads) != 0:
                os.system(bedtools + ' intersect -a %s -b %s -wa | sort -u > %s' % (reads, v2nBed_break, newReads))
                if os.path.exists(newReads) and os.path.getsize(newReads) != 0:
                    readsDF = pd.read_csv(newReads, sep='\t', names=['ins', 'start', 'end', 'read', 'MAPQ', 'strand'])
            if os.path.exists(split) and os.path.getsize(split) != 0:
                splitDF = pd.read_csv(split, sep='\t', names=['ins', 'start', 'end', 'read'])
            INSid = []
            if 'splitDF' in locals().keys():
                INSid.extend(list(np.unique(list(splitDF['ins']))))
            if 'readsDF' in locals().keys():
                INSid.extend(list(np.unique(list(readsDF['ins']))))
            INSid = list(set(INSid))
            for ins in INSid:
                ins_info = ins.split("_")
                if 'readsDF' in locals().keys():
                    sub_reads = readsDF[readsDF['ins'] == ins]
                    if not sub_reads.empty:
                        ins_info.append(str(len(np.unique(sub_reads['read']))))
                    else:
                        ins_info.append('0')
                else:
                    ins_info.append('0')
                if 'splitDF' in locals().keys():
                    sub_split = splitDF[splitDF['ins'] == ins]
                    if not sub_split.empty:
                        ins_info.append(str(len(np.unique(sub_split['read']))))
                    else:
                        ins_info.append('0')
                else:
                    ins_info.append('0')
                newline = '\t'.join(ins_info) + '\n'
                reads_result.write(newline)
            if 'splitDF' in locals().keys():
                del splitDF
            if 'readsDF' in locals().keys():
                del readsDF
    reads_result.close()
    out = open('insExtend_site_samples_reads.%s_INS_stat.bed' % idx, 'w')  # 要求nonsplit read必须在INS之内或者跨过
    for line in open("samples_reads_INS_stat_%s.txt" % idx).readlines():
        if 'start' not in line:
            info = line.strip().split('\t')
            loc_refname = '_'.join(info[0:5])
            if ('C1' in info[0]) or ('C2' in info[0]) or ('C3' in info[0]):
                sample = info[0]
                Sample = sample
            else:
                sample = info[0].split("C")[0]
                Sample = sample + 'C'
            refname = Sample + '_' + '_'.join(info[1:5])
            hg38_loc = open(
                os.path.join(ins_locdir, Sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')).readline()
            if '-1' not in hg38_loc:
                newline = hg38_loc.strip() + '\t' + '\t'.join(info[0:7]) + '\n'
                out.write(newline)
    out.close()
    os.system('sort -k1,1 -k2,2n insExtend_site_samples_reads.%s_INS_stat.bed > \
                sorted_samples_reads.%s_INS_stat.bed' % (idx, idx))
    os.system(bedtools + ' closest -a sorted_samples_reads.%s_INS_stat.bed \
                -b %s -t all -D ref > gene_samples_reads.%s_INS_stat.bed' % (idx, genebed, idx))
    os.system(bedtools + ' intersect -a sorted_samples_reads.%s_INS_stat.bed \
                -b %s -wa -wb | sort -k1,1 -k2,2n -u > \
                gene5kb_overlap_samples_reads.%s_INS_stat.bed' % (idx, gene5kb_bed, idx))

## 确定split segment 比对位置在ref的位置
def new_segment_gene(sample, idx):
    print(sample)
    if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
        Sample = sample
    else:
        Sample = sample + 'C'
    inssplit = os.path.join(workdir, sample,
                            sample + '_ins_split-reads.%s_detail.txt' % idx)
    allsplit = os.path.join(os.path.join(workdir, sample, 'sorted_uniq_' + sample +
                                         '_ins_split-reads.%s_detail.bed' % idx))
    if os.path.getsize(allsplit) > 0 and os.path.getsize(inssplit) > 0:
        INSsplit_bed = pd.read_csv(inssplit,
                                   sep='\t', header=None, names=['refname', 'start', 'end', 'read'])
        ALLsplit_bed = pd.read_csv(allsplit,
                                   sep='\t', header=None, names=['refname', 'start', 'end', 'read'])
        split_bed = ALLsplit_bed[ALLsplit_bed['read'].isin(list(INSsplit_bed['read']))]
        new_split_bed = open(os.path.join(workdir, sample,
                                          'hg38_' + sample + '_ins_split-reads.%s_detail.bed' % idx), 'w')
        split_bed[['start', 'end']] = split_bed[['start', 'end']].astype(str)
        for refname in np.unique(list(split_bed['refname'])):
            loc_hg38 = os.path.join(ins_locdir, Sample, refname,
                                    refname + '_SRM_INS_hg38_coordinate.bed')
            v2n_loc_file = os.path.join(ins_locdir, Sample, refname, refname + '_SRM_INS_coordinate.bed')
            if not os.path.exists(loc_hg38):
                loc_hg38 = os.path.join(dup_locdir, Sample, refname,
                                        refname + '_SRM_INS_hg38_coordinate.bed')
                v2n_loc_file = os.path.join(dup_locdir, Sample, refname, refname + '_SRM_INS_coordinate.bed')
            hg38_loc = pd.read_csv(loc_hg38, sep='\t', header=None)
            hg38Start = hg38_loc.iloc[0, 1]
            hg38End = hg38_loc.iloc[0, 2]
            hg38Chr = hg38_loc.iloc[0, 0]
            if (str(hg38End) != '-1') and (str(hg38Start) != '-1'):
                v2n_loc = pd.read_csv(v2n_loc_file, sep='\t', header=None)
                v2nStart = v2n_loc.iloc[0, 1]
                v2nEnd = v2n_loc.iloc[0, 2]
                sub_split_bed = split_bed[split_bed['refname'] == refname]
                sub_split_bed[['start', 'end']] = sub_split_bed[['start', 'end']].astype(str)
                for i in range(0, sub_split_bed.shape[0]):
                    if int(sub_split_bed.iloc[i, 2]) < v2n_loc.iloc[0, 1]:
                        newstart = hg38Start - v2nStart + int(sub_split_bed.iloc[i, 1])
                        newend = hg38Start - v2nStart + int(sub_split_bed.iloc[i, 2])
                        newline = '\t'.join([hg38Chr, str(newstart), str(newend)]) + '\t' + '\t'.join(
                            list(sub_split_bed.iloc[i,])) + '\n'
                    elif int(sub_split_bed.iloc[i, 1]) > v2n_loc.iloc[0, 2]:
                        newstart = hg38End - v2nEnd + int(sub_split_bed.iloc[i, 1])
                        newend = hg38End - v2nEnd + int(sub_split_bed.iloc[i, 2])
                        newline = '\t'.join([hg38Chr, str(newstart), str(newend)]) + '\t' + '\t'.join(
                            list(sub_split_bed.iloc[i,])) + '\n'
                    else:
                        sites = [int(sub_split_bed.iloc[i, 1]), int(sub_split_bed.iloc[i, 2]),
                                 v2n_loc.iloc[0, 1], v2n_loc.iloc[0, 2]]
                        sites.sort()
                        if sites == [int(sub_split_bed.iloc[i, 1]), v2n_loc.iloc[0, 1],
                                     int(sub_split_bed.iloc[i, 2]), v2n_loc.iloc[0, 2]]:
                            newstart = hg38Start - v2nStart + int(sub_split_bed.iloc[i, 1])
                            newend = hg38Start
                            newline = '\t'.join([hg38Chr, str(newstart), str(newend)]) + '\t' + '\t'.join(
                                list(sub_split_bed.iloc[i,])) + '\n'
                        elif sites == [v2n_loc.iloc[0, 1], int(sub_split_bed.iloc[i, 1]),
                                       v2n_loc.iloc[0, 2], int(sub_split_bed.iloc[i, 2])]:
                            newstart = hg38End
                            newend = hg38End - v2nEnd + int(sub_split_bed.iloc[i, 2])
                            newline = '\t'.join([hg38Chr, str(newstart), str(newend)]) + '\t' + '\t'.join(
                                list(sub_split_bed.iloc[i,])) + '\n'
                        else:  # segment正好落在INS区域内，hg38上没有对应位置
                            newline = ''
                    if any(newline):
                        new_split_bed.write(newline)
        new_split_bed.close()
        gene_new_split_bed = os.path.join(workdir, sample, "gene_hg38_" + sample +
                                          '_ins_split-reads.%s_detail.bed' % idx)
        New_split_bed = os.path.join(workdir, sample, 'hg38_' + sample +
                                     '_ins_split-reads.%s_detail.bed' % idx)
        exon_new_split_bed = os.path.join(workdir, sample, "exon_hg38_" + sample +
                                          '_ins_split-reads.%s_detail.bed' % idx)
        os.system(bedtools + ' intersect -a %s -b %s -wa -wb > %s' % (New_split_bed, genebed, gene_new_split_bed))
        os.system(bedtools + ' intersect -a %s -b %s -wa -wb > %s' % (New_split_bed, exonbed, exon_new_split_bed))


pools = Pool(5)
for sample in list(sample_list.keys()):
    for idx in ['all-split', 'disc']:
        pools.apply_async(new_segment_gene, args=(sample, idx,))

pools.close()
pools.join()
del pools

####### merge split/disc/unmapped reads
workdir = '/nanopore/xl/INS_DUP_RNA'
ins_locdir = '/Assembly_result'
dup_locdir = '/Assembly_DUP/DUP_2286/Assembly_result'
os.chdir(workdir)
oldresult = pd.read_csv(os.path.join(workdir, 'allSamples_repeat_gene_INS_stat.txt'), sep='\t')
newresult = oldresult.copy()
for idx in ['disc', 'all-split']:
    reads = pd.read_csv('samples_reads_INS_stat_%s.txt' % idx, sep='\t')
    reads['newid_xl'] = reads['sample'] + '_' + reads['chr'] + '_' + reads['start'].map(str) \
                        + '_' + reads['end'].map(str) + '_' + reads['id'].map(str)
    reads.columns = ['sample', 'chr', 'start', 'end', 'id', 'disc-reads', 'disc-split_reads', 'newid_xl']
    gene5kb = pd.read_csv("gene5kb_overlap_samples_reads.%s_INS_stat.bed" % idx, sep='\t', header=None)
    gene5kb1 = gene5kb[[3, 4, 5, 6, 7, 8, 9, 15, 16]]
    gene5kb1.columns = ['sample', 'chr', 'start', 'end', 'id', idx + '_read', idx + '_split-read', 'genetype', 'gene']
    gene5kb1['newid_xl'] = gene5kb1['sample'] + '_' + gene5kb1['chr'] + '_' + gene5kb1['start'].map(str) + \
                           '_' + gene5kb1['end'].map(str) + '_' + gene5kb1['id'].map(str)
    gene5kb2 = gene5kb1[['gene', 'genetype', 'newid_xl', idx + '_read', idx + '_split-read']]
    os.system('cat */gene_hg38_Lung*_ins_split-reads.%s_detail.bed > \
                    allSamples_gene_hg38_split-reads.%s_detail.txt' % (idx, idx))
    all_repeat_gene = pd.read_csv("allSamples_gene_hg38_split-reads.%s_detail.txt" % idx, header=None,
                                  sep='\t')
    all_repeat_gene1 = all_repeat_gene[[3, 6, 12, 13]]
    all_repeat_gene1.columns = ['svid', 'readid', 'genetype', 'gene']
    all_repeat_gene1['sample'] = all_repeat_gene1['svid'].str.split("_").str[0]
    all_repeat_gene2 = all_repeat_gene1.drop_duplicates()
    all_repeat_gene3 = all_repeat_gene2.pivot_table(index=['svid', 'gene', 'genetype', 'sample'],
                                                    values='readid',
                                                    aggfunc={'readid': 'count'},
                                                    fill_value=0)
    all_repeat_gene3['gene'] = all_repeat_gene3.index.get_level_values('gene')
    all_repeat_gene3['sample'] = all_repeat_gene3.index.get_level_values('sample')
    all_repeat_gene3['newid_xl'] = all_repeat_gene3.index.get_level_values('svid')
    all_repeat_gene3['genetype'] = all_repeat_gene3.index.get_level_values('genetype')
    all_repeat_gene3 = all_repeat_gene3[['gene', 'genetype', 'newid_xl', 'readid']]
    all_repeat_gene3.columns = ['gene', 'genetype', 'newid_xl', idx + '_split-read-gene']
    all_repeat_gene3.index = list(range(0, all_repeat_gene3.shape[0]))
    newdf = pd.merge(gene5kb2, all_repeat_gene3, on=['gene', 'genetype', 'newid_xl'], how='outer')
    newdf1 = newdf.fillna(0)
    newdf1['sample'] = newdf1['newid_xl'].str.split("_").str[0]
    newresult = pd.merge(newresult, newdf1, on=['gene', 'genetype', 'newid_xl', 'sample'], how='outer')

dup_type_all = pd.read_csv(workdir + "/DUP_repeat.bed", header=None, sep='\t')
ins_type_all = pd.read_csv(workdir + "/INS_repeat.bed", header=None, sep='\t')

newresult1 = newresult.fillna(0)
for i in range(0, newresult1.shape[0]):
    if newresult1.iloc[i, 1] == 0:
        Sample = newresult1.iloc[i, 5].split("_")[0]
        refname = newresult1.iloc[i, 5]
        hg38_loc_path = os.path.join(ins_locdir, Sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
        if not os.path.exists(hg38_loc_path):
            hg38_loc_path = os.path.join(dup_locdir, Sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
        hg38_loc = open(hg38_loc_path).readline()
        if '-1' not in hg38_loc:
            newresult1.iloc[i, 0], newresult1.iloc[i, 1], newresult1.iloc[i, 2] = hg38_loc.strip().split()
        if refname in list(dup_type_all[0]):
            newresult1.iloc[i, 10], newresult1.iloc[i, 11] = list(
                dup_type_all[dup_type_all[0] == refname].iloc[0, [3, 4]])
        elif refname in list(ins_type_all[0]):
            newresult1.iloc[i, 10], newresult1.iloc[i, 11] = list(
                ins_type_all[ins_type_all[0] == refname].iloc[0, [3, 4]])
    if newresult1.iloc[i, 12] == 0:
        teList = []
        for info in newresult1.iloc[i, 10].split("|"):
            if 'mask' not in info:
                teList.append('_'.join(info.split("_")[2:len(info.split("_"))]))
            else:
                teList.append(info)
        newresult1.iloc[i, 12] = ','.join(list(set(teList)))

newresult1.columns = ['#chr', 'start', 'end', 'gene', 'genetype', 'newid_xl', 'unmapped_nonsplit-read',
                      'unmapped_split-read', 'unmapped_split_read_gene', 'sample', 'repeat', 'repeat_class',
                      'repeatType', 'disc_nonsplit-read', 'disc_split-read', 'disc_split-read-gene',
                      'all-split_nonsplit-read', 'all-split_split-read', 'all-split_split-read-gene']
newresult1 = newresult1[['#chr', 'start', 'end', 'gene', 'genetype', 'newid_xl', 'sample', 'repeat', 'repeat_class',
                         'repeatType', 'unmapped_nonsplit-read', 'unmapped_split-read', 'unmapped_split_read_gene',
                         'disc_nonsplit-read', 'disc_split-read', 'disc_split-read-gene',
                         'all-split_nonsplit-read', 'all-split_split-read', 'all-split_split-read-gene']]
newresult1['start'] = newresult1['start'].map(int)
newresult1['end'] = newresult1['end'].map(int)
newresult1.to_csv("allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", sep='\t', index=0)
teinfo2 = teinfo1[['newid_xl', 'segment_dup', 'repeat_dup', 'total_len']]
newresult2 = pd.merge(newresult1, teinfo2, on=['newid_xl'], how='left')
ins_newresult = newresult2[(newresult2['segment_dup'] == 0) & (newresult2['repeat_dup'] == 0)]
dup_newresult = newresult2[(newresult2['segment_dup'] != 0) | (newresult2['repeat_dup'] != 0)]
ins_newresult.to_csv("INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", sep='\t', index=0)
dup_newresult.to_csv("TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", sep='\t', index=0)

anno = {'gene': '/ref/gencode.v33.annotation_gene_20200312.bed',
        'exon': '/ref/gencode.v33.annotation_exon_20200312.bed',
        'promoter': '/ref/gencode.v33.annotation_promoter_20200312.bed',
        'intron': '/ref/gencode.v33.annotation_intron_20200312.bed',
        'downstream': '/ref/gencode.v33.annotation_gene_downstream2kb.bed',
        'intergenic': '/ref/gencode.v33.annotation_intergenic_20200312.bed',
        'enhancer': ' /ref/human_permissive_enhancers_phase_1_and_2_LIFTOVER_hg38_new.bed',
        'openPeak': '/data/fs08 /Assembly_ztf/ATAC_peak/inter_samples_summits_merged_filtered.bed'}

file_list = ['allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt',
             "INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt",
             "TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt"]
for file in file_list:
    for dis in [100, 200, 500]:
        all = pd.read_csv(allINS_rna, sep='\t')
        newbed = 'result/%s' % file.split('')[0] + str(dis) + '.bed'
        all['start'] = all['start'] - dis
        all['end'] = all['end'] + dis
        all.to_csv(newbed, sep='\t', index=0, header=0)

total_ins_rna = '1775'
totalINS = '1423'
totalTD = '352'
new_file_list = {'all': ['allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt', total_ins_rna],
                 'INS': ["INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", totalINS],
                 'TD': ["TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", totalTD]}
outfile = open('allSamples_repeat_gene_INS_stat_unmapped_disc_split_for-Oddratio.txt', 'w')
for annoregion in anno:
    for idx in new_file_list:
        os.system(bedtools + ' intersect -a %s -b %s -wa -wb > result/%s_%s' %
                  (new_file_list[idx][0], anno[annoregion], annoregion, new_file_list[idx][0]))
        count = os.popen("awk -F '\t' '{print $6}' result/%s_%s \
                            | sort -u | wc -l " % (annoregion, allINS_rna)).read().strip().split()[0]
        newline = '\t'.join([idx, annoregion, '0', new_file_list[idx][1], count]) + '\n'
        outfile.write(newline)
        for dis in [100, 200, 500]:
            anno_result = 'result/%s' % new_file_list[idx][0].split('')[0] + \
                          str(dis) + '.bed'
            os.system(bedtools + ' intersect -a %s -b %s -wa -wb > result/%s_%s' % (
                anno_result, anno[annoregion], annoregion, anno_result.split("/")[-1]))
            count = os.popen("awk -F '\t' '{print $6}' result/%s_%s \
                                    | sort -u | wc -l " % (
            annoregion, anno_result.split("/")[-1])).read().strip().split()[0]
            newline = '\t'.join([idx, annoregion, str(dis), new_file_list[idx][1], count]) + '\n'
            outfile.write(newline)

outfile.close()



#######################################################################################################################
##########
# all RNA reads mapped to candidate INS sequence
gene_list = os.popen("grep -v '#' INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt | \
                        awk -F '\t' '{print $4}' | sort -u").read().strip().split('\n')
insbed = ' /new_somatic_SV_TGS_bed/somatic_sv_new_bed/all_samples_newINS_candidateRegion.bed'
os.system(bedtools + ' intersect -a %s -b %s -wa -wb > gene5kb_all_samples_newINS_candidateRegion.txt' % (
gene5kb_bed, insbed))
gene5kb_ins = pd.read_csv("gene5kb_all_samples_newINS_candidateRegion.txt", sep='\t', header=None)
gene5kb_ins1 = gene5kb_ins[gene5kb_ins[6].isin(gene_list)]
gene5kb_ins1['ins_region'] = gene5kb_ins1[11] + '_' + gene5kb_ins1[12].map(str) + '_' + gene5kb_ins1[13].map(str) +\
                             '_' + gene5kb_ins1[15]
teinfo1['sample'] = teinfo1['newid_xl'].str.split("_").str[0]
teinfo1['new_ins_region'] = teinfo1['ins_region'] + '_' + teinfo1['sample']
ins_teinfo1_gene5kb = teinfo1[teinfo1['new_ins_region'].isin(list(gene5kb_ins1['ins_region']))]
rna_fq_dir = '/data/fs02 /nanopore_wxy/RNA-seq/cleandata'

resultdir = os.path.join(workdir, 'assembled_candidate')
if not os.path.exists(resultdir):
    os.mkdir(resultdir)


def mapping_assemble(sample, refname):
    sub_resultdir = os.path.join(workdir, 'assembled_candidate', sample)
    if not os.path.exists(sub_resultdir):
        os.mkdir(sub_resultdir)
    os.chdir(sub_resultdir)
    if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
        Sample = sample
    else:
        Sample = sample.split("C")[0]
    fq1 = os.path.join(rna_fq_dir, Sample, Sample, Sample + '_1.fq.gz')
    fq2 = os.path.join(rna_fq_dir, Sample, Sample, Sample + '_2.fq.gz')
    ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
    locfile = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
    if not os.path.exists(ins_extend):
        ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
        locfile = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
    os.system('cp %s ./' % locfile)
    os.system(hisat2_build + " %s %s -p 10" % (ins_extend, sub_resultdir))
    sam = sample + '_allRNAreads_mapped_%s.sam' % refname
    os.system(hisat2 + ' -p 10 -x {sample_workdir} -1 {fq1} -2 {fq2} -S {sam}'.format(
        sample_workdir=sub_resultdir, fq1=fq1, fq2=fq2, sam=sam))
    os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (sam, sam.split(".sam")[0]))
    os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (sam.split(".sam")[0], sam.split(".sam")[0]))
    os.system('rm *sam')
    os.system(stringtie + " sorted_%s.bam -m 100 -c 1 -o %s.gtf -p 10" % (sam.split(".sam")[0], sam.split(".sam")[0]))


pools = Pool(5)
for sample in list(set(ins_teinfo1_gene5kb['sample_id'])):
    subdf = ins_teinfo1_gene5kb[ins_teinfo1_gene5kb['sample_id'] == sample]
    for i in range(0, subdf.shape[0]):
        refname = subdf.iloc[i, 0]
        pools.apply_async(mapping_assemble, args=(sample, refname,))

pools.close()
pools.join()
del pools



##########
# gene structural 对应到 assembled INS上
ins_rna = pd.read_csv("/INS_DUP_RNA/INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split.txt", sep='\t')
# os.system('grep -v "#"   /ref/gencode.v33.annotation.gtf >   /ref/gencode.v33.annotation_noanno.gtf')
gene_bed = pd.read_csv('/ref/gencode.v33.annotation_noanno.gtf', sep='\t', header=None)
gene_bed['gene'] = gene_bed[8].str.split('gene_name "').str[1]
gene_bed['gene'] = gene_bed['gene'].str.split('"').str[0]
gene_bed['gene_id'] = gene_bed[8].str.split('gene_id "').str[1]
gene_bed['gene_id'] = gene_bed['gene_id'].str.split('"').str[0]

ins_rna1 = ins_rna[['gene', 'newid_xl', 'total_len']]
ins_rna1_gene = pd.merge(ins_rna1, gene_bed, on=['gene'], how='left')
newgtf = open('INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf', 'w')
exon_newgtf = open('exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf', 'w')
if os.path.exists("newgtf.log"):
    os.system('rm newgtf.log')


for i in range(0, ins_rna1_gene.shape[0]):
    info = ins_rna1_gene.iloc[i].map(str).tolist()
    refname = info[1]
    ins_hg38_bam = os.path.join(ins_locdir, refname.split('_')[0], refname, refname + '_SRM2hg38.bam')
    ins_loci = os.path.join(ins_locdir, refname.split('_')[0], refname, refname + '_SRM_INS_coordinate.bed')
    if not os.path.exists(ins_hg38_bam):
        ins_hg38_bam = os.path.join(dup_locdir, refname.split('_')[0], refname, refname + '_SRM2hg38.bam')
        ins_loci = os.path.join(ins_locdir, refname.split('_')[0], refname, refname + '_SRM_INS_coordinate.bed')
    ins_loci_df = pd.read_csv(ins_loci, sep='\t',header=None)
    bamfile = pysam.AlignmentFile(ins_hg38_bam, "rb")
    allreads = bamfile.fetch()
    reads1 = next(allreads)
    align_info = reads1.get_aligned_pairs(matches_only=False, with_seq=False)
    del_idx = []
    for m in range(0, len(align_info)):
        if align_info[m][1] is None:
            del_idx.append(m)
        if align_info[m][0] is None:
            del_idx.append(m)
    if del_idx:
        counter = 0
        for index in del_idx:
            index = index - counter
            align_info.pop(index)
            counter += 1
    bamfile.close()
    newinfo = info[3:len(info)]
    newinfo[0] = refname
    if (align_info[0][1] < int(newinfo[3])) & (align_info[-1][1] < int(newinfo[3])):
        next_step = False
    elif (align_info[0][1] > int(newinfo[4])) & (align_info[-1][1] > int(newinfo[4])):
        next_step = False
    else:
        next_step = True
    if next_step:
        for j in range(3, 5):
            site = int(newinfo[j])
            for n in range(0, len(align_info)):
                pair = align_info[n]
                if site - pair[1] in range(0, 21):
                    print(pair)
                    newinfo[j] = str(pair[0]+1)
                    idx = True
            if "idx" not in locals():
                idx = False
            if not idx:
                if abs(align_info[0][1] - site) < abs(align_info[-1][1] - site):
                    newinfo[j] = str(align_info[0][0]+1)
                else:
                    newinfo[j] = str(align_info[-1][0]+1)
                newinfo[-1] = newinfo[-1] + 'type "incomplete";'
            else:
                newinfo[-1] = newinfo[-1] + 'type "complete";'
            del idx
        if int(newinfo[3]) <= int(newinfo[4]):
            newline = '\t'.join(newinfo) + '\n'
            newgtf.write(newline)
            if newinfo[2]=='exon':
                exon_newgtf.write(newline)
        else:
            os.system('echo "%s" >> newgtf.log' % ('\t'.join(info)))
            print("error")
            break

newgtf.close()
exon_newgtf.close()
os.system('sort -k1,1 -k4,4n -k5,5n INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf > \
            sorted_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf')
os.system('sort -k1,1 -k4,4n -k5,5n exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf > \
            sorted_exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_gencodeV33.gtf')

############
te_bed = os.path.join(workdir, 'all.INS.DUP.anno.region.clean.te_list.bed')
genebed = '  /ref/sorted_gencode.v33.annotation_gene_20200312.bed'

gene_new_split_bed_disc = 'gene_hg38_allSample_ins_split-reads.dics_detail.bed'
os.system('cat '+workdir+"/*/gene_hg38_*_ins_split-reads.disc_detail.bed | sort -k1,1 -k2,2n > " +
          gene_new_split_bed_disc)
gene_new_split_bed_all_split = 'gene_hg38_allSample_ins_split-reads.all-split_detail.bed'
os.system('cat '+workdir+"/*/gene_hg38_*_ins_split-reads.all-split_detail.bed | sort -k1,1 -k2,2n > " +
          gene_new_split_bed_all_split)
gene_new_split_bed_unmapped = 'gene_hg38_allSample_ins_split-reads.unmapped_detail.bed'
os.system('cat '+workdir+'/*/gene_hg38_*_ins_split-reads_detail_distance.bed| sort -k1,1 -k2,2n > ' +
          gene_new_split_bed_unmapped)
gene_hg38_reads = {'disc':gene_new_split_bed_disc,
                   'all-split':gene_new_split_bed_all_split,
                   'unmapped':gene_new_split_bed_unmapped}

unmapped_newReads = os.path.join(workdir, 'allSample_ins_reads_nonsplit_in_INS.txt')
os.system('cat %s/*/*_ins_reads_nonsplit_in_INS.txt | sort -k1,1 -k2,2n > %s' % (workdir, unmapped_newReads))
unmapped_split = os.path.join(workdir, 'allSample_ins_split-reads_detail.txt')
os.system('cat %s/*/*_ins_split-reads_detail.txt | sort -k1,1 -k2,2n > %s' % (workdir, unmapped_split))
for idx in ['all-split', 'disc']:
    newReads = os.path.join(workdir, 'allSample_ins_reads.%s_nonsplit_in_INS.txt' % idx)
    os.system('cat %s/*/*_ins_reads.%s_nonsplit_in_INS.txt | sort -k1,1 -k2,2n > %s' % (
                workdir, idx, newReads))
    split = os.path.join(workdir, 'allSample_ins_split-reads.%s_detail.txt' % idx)
    os.system('cat %s/*/*_ins_split-reads.%s_detail.txt | sort -k1,1 -k2,2n > %s' % (
                workdir, idx, split))

reads_list = {'unmapped_nonsplit-read':'allSample_ins_reads_nonsplit_in_INS.txt',
              'unmapped_split-read':'allSample_ins_split-reads_detail.txt',
              'disc_nonsplit-read':'allSample_ins_reads.disc_nonsplit_in_INS.txt',
              'disc_split-read':'allSample_ins_split-reads.disc_detail.txt',
              'all-split_nonsplit-read':'allSample_ins_reads.all-split_nonsplit_in_INS.txt',
              'all-split_split-read':'allSample_ins_split-reads.all-split_detail.txt'}
for read in list(reads_list.keys()):
    os.system(bedtools + ' intersect -a %s -b %s -wa -wb | grep "\<INS\>" > te_INS_%s' % (reads_list[read],te_bed, reads_list[read]))
    te_result = pd.read_csv('te_INS_' + reads_list[read],sep='\t',header=None)
    if read in ['disc_split-read','unmapped_split-read','all-split_split-read']:
        te_result1 = te_result[[0, 3, 7,8,9, 10, 11]]
    else:
        te_result1 = te_result[[0,3,9,10,11,12,13]]
    te_result1.columns = ['refname','read','SV','sample','te_type','te','te_class']
    te_result2 = pd.pivot_table(te_result1, index=['refname','sample','te_type','te','te_class'],
                                values=['read'], aggfunc=[len])
    te_result3 = te_result2.reset_index()
    te_result3.columns = ['refname','sample','te_type','te','te_class','read_count']
    te_result3['read_type'] = [read]*te_result3.shape[0]
    if '_split-read' in read:
        idx = read.split("_")[0]
        genebed = gene_hg38_reads[idx]
        genebed_df = pd.read_csv(genebed,sep='\t',header=None)
        genebed_df1 = genebed_df[genebed_df[3].isin(list(te_result3['refname']))]
        gene_list = []
        gene_type_list = []
        for i in range(0, te_result3.shape[0]):
            if list(te_result3['refname'])[i] in list(genebed_df1[3]):
                gene_list.append(';'.join(list(set(genebed_df1[genebed_df1[3]==list(te_result3['refname'])[i]][13]))))
                gene_type_list.append(';'.join(list(set(genebed_df1[genebed_df1[3] == list(te_result3['refname'])[i]][12]))))
            else:
                gene_list.append('')
                gene_type_list.append('')
        te_result3['gene'] = gene_list
        te_result3['gene_type'] = gene_type_list
    else:
        gene_list = []
        gene_type_list = []
        for i in range(0, te_result3.shape[0]):
            if list(te_result3['refname'])[i] in list(ins_rna['newid_xl']):
                gene_list.append(';'.join(list(set(ins_rna[ins_rna['newid_xl'] ==
                                                           list(te_result3['refname'])[i]]['gene']))))
                gene_type_list.append(
                    ';'.join(list(set(ins_rna[ins_rna['newid_xl'] ==
                                              list(te_result3['refname'])[i]]['genetype']))))
            else:
                gene_list.append('')
                gene_type_list.append('')
        te_result3['gene'] = gene_list
        te_result3['gene_type'] = gene_type_list
    te_result3.to_csv('te_stat/stat_te_INS_' + reads_list[read], sep='\t', index=0)

