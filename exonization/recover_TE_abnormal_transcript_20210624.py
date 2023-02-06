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

workdir = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524//INS_DUP_RNA_20210624'
if not os.path.exists(workdir):
    os.mkdir(workdir)

os.chdir(workdir)
hg38_resultDir = '/data/fs08/wangzf/nanopore/rna/RNA-pipeline_bam'
ins_locdir = '/data/fs08/wangzf/Assembly_ztf/INS_RepeatMasker_f50bp/Assembly_result'
dup_locdir = '/data/fs08/wangzf/Assembly_DUP/DUP_2286/Assembly_result'

bedtools = '/data/fs01/biosoft/bedtools-2.28.0/bin/bedtools'
samtools = '/data/fs01/biosoft/samtools-1.9/samtools'
hisat2 = '/data/fs01/wangzf/software/hisat2-2.1.0/hisat2'
hisat2_build = '/data/fs01/wangzf/software/hisat2-2.1.0/hisat2-build'
stringtie = '/data/fs01/wangzf/software/stringtie-2.1.4/stringtie'
stringtie_2_1_6 = '/data/fs08/wangzf/nanopore/xl/tools/stringtie-2.1.6.Linux_x86_64/stringtie'
######################################################################################################
## 制作TE bed
te_info = '/data/fs08/wangzf/nanopore/pxn/INS/all.INS.DUP.anno.region.clean.20210622.tsv'
teinfo = pd.read_csv(te_info, sep='\t')
teinfo1 = teinfo[['newid_xl', 'total_len', 'clean_decision', 'sample_id', 'ins_region', 'decision_type',
                  'segment_dup', 'repeat_dup']]
dup_info = teinfo1[(teinfo1['segment_dup'] == 1) | (teinfo1['repeat_dup'] == 1)]  # 5519
ins_info = teinfo1[~teinfo1.index.isin(list(dup_info.index))]  # 5640
# dup
dup_type = open("DUP_repeatType_20210624.bed", 'w')
dup_type_all = open("DUP_repeat_20210624.bed", 'w')
for i in range(0, dup_info.shape[0]):
    dupinfo = dup_info.iloc[i].values.tolist()
    sample = dupinfo[3]
    locfile = os.path.join(ins_locdir, dupinfo[3], dupinfo[0], dupinfo[0] + '_SRM_INS_coordinate.bed')
    if os.path.exists(locfile):
        loc = pd.read_csv(locfile, sep='\t', header=None)
    else:
        locfile = os.path.join(dup_locdir, dupinfo[3], dupinfo[0], dupinfo[0] + '_SRM_INS_coordinate.bed')
        if os.path.exists(locfile):
            loc = pd.read_csv(locfile, sep='\t', header=None)
        else:
            print('error ' + dupinfo[0])
            break
    type_list = dupinfo[2].split("|")
    newloc = loc.copy()
    newloc.iloc[0, 1] = newloc.iloc[0, 1] - 50
    newloc.iloc[0, 2] = newloc.iloc[0, 2] + 50
    for type in type_list:
        newinfo = newloc.iloc[0].map(str).values.tolist()
        newinfo[0] = dupinfo[0]
        if type != 'Unmask':
            type_info = type.split("_")
            start = int(type_info[0]) - 1
            end = int(type_info[1]) - 1
            Type = '_'.join(type_info[2:len(type_info)])
            newinfo[1] = str(int(newinfo[1]) + start)
            newinfo[2] = str(int(newinfo[1]) + end)
            newinfo[-1] = Type
        else:
            newinfo[-1] = type
        newline = '\t'.join(newinfo) + '\n'
        dup_type.write(newline)
    newloc.iloc[0, 0] = dupinfo[0]
    newloc.iloc[0, 3] = dupinfo[2]
    if len(type_list) > 1:
        new_line = '\t'.join(newloc.iloc[0].map(str).tolist()) + '\tmultiple\n'
    else:
        new_line = '\t'.join(newloc.iloc[0].map(str).tolist()) + '\tsingle\n'
    dup_type_all.write(new_line)

dup_type_all.close()
dup_type.close()
# ins
ins_type = open("INS_repeatType_20210624.bed", 'w')
ins_type_all = open("INS_repeat_20210624.bed", 'w')
for i in range(0, ins_info.shape[0]):
    insinfo = ins_info.iloc[i].values.tolist()
    sample = insinfo[3]
    locfile = os.path.join(ins_locdir, insinfo[3], insinfo[0], insinfo[0] + '_SRM_INS_coordinate.bed')
    if os.path.exists(locfile):
        loc = pd.read_csv(locfile, sep='\t', header=None)
    else:
        locfile = os.path.join(dup_locdir, insinfo[3], insinfo[0], insinfo[0] + '_SRM_INS_coordinate.bed')
        if os.path.exists(locfile):
            loc = pd.read_csv(locfile, sep='\t', header=None)
        else:
            print('error ' + insinfo[0])
            break
    type_list = insinfo[2].split("|")
    newloc = loc.copy()
    newloc.iloc[0, 1] = newloc.iloc[0, 1] - 50
    newloc.iloc[0, 2] = newloc.iloc[0, 2] + 50
    for type in type_list:
        newinfo = newloc.iloc[0].map(str).values.tolist()
        newinfo[0] = insinfo[0]
        if type != 'Unmask':
            type_info = type.split("_")
            start = int(type_info[0]) - 1
            end = int(type_info[1]) - 1
            Type = '_'.join(type_info[2:len(type_info)])
            newinfo[1] = str(int(newinfo[1]) + start)
            newinfo[2] = str(int(newinfo[1]) + end)
            newinfo[-1] = Type
        else:
            newinfo[-1] = type
        newline = '\t'.join(newinfo) + '\n'
        ins_type.write(newline)
    newloc.iloc[0, 0] = insinfo[0]
    newloc.iloc[0, 3] = insinfo[2]
    if len(type_list) > 1:
        new_line = '\t'.join(newloc.iloc[0].map(str).tolist()) + '\tmultiple\n'
    else:
        new_line = '\t'.join(newloc.iloc[0].map(str).tolist()) + '\tsingle\n'
    ins_type_all.write(new_line)

ins_type.close()
ins_type_all.close()

# the site of TE on INS-assembled sequence (bed)
tf_te_info = pd.read_csv("/data/fs08/wangzf/nanopore/ztf/INSandDUP_V3/INS_11159/11159_INS_reannotate_duplication_repeatinfo.tsv",
                         sep='\t',header=None)
tf_te_info['type'] = ['DUP']*tf_te_info.shape[0]
tf_te_info1 = pd.read_csv("/data/fs08/wangzf/nanopore/ztf/INSandDUP_V3/INS_11159/11159_INS_reannotate_insertion_repeatinfo.tsv",
                         sep='\t',header=None)
tf_te_info1['type'] = ['INS']*tf_te_info1.shape[0]
tf_te_info = tf_te_info.append(tf_te_info1)
out = open("all.INS.DUP.anno.region.clean.20210622.te_list.bed",'w')
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
        else: # Lung100C_chr2_3170991_3171145_41-INS
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
for Sample in os.listdir("/data/fs08/wangzf/nanopore/rna/extract_hg38/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/data/fs08/wangzf/nanopore/rna/extract_hg38/' + Sample
            sample_list[Sample] = sampledir

for Sample in os.listdir("/data/fs08/wangzf/nanopore/rna/extract_hg38_1/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/data/fs08/wangzf/nanopore/rna/extract_hg38_1/' + Sample
            sample_list[Sample] = sampledir

for Sample in os.listdir("/data/fs08/wangzf/nanopore/rna/extract_hg38_2/"):
    if 'ht2' not in Sample:
        if 'Lung' in Sample:
            sampledir = '/data/fs08/wangzf/nanopore/rna/extract_hg38_2/' + Sample
            sample_list[Sample] = sampledir


# 提取split reads
def split_reads(sample, bam):
    log = open(sample + '_split.log', 'w')
    log.write(sample + " " + time.asctime() + ' start\n')
    os.system(samtools + ' index ' + bam)
    samfile = pysam.AlignmentFile(bam, 'rb')
    outbam = pysam.AlignmentFile(sample + '_split.bam', 'wb', template=samfile)
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
    split_reads.close()
    log.write(sample + " " + time.asctime() + ' end\n')
    log.close()
    os.system('sort -u %s_split_readsname.txt > sorted_uniq_%s_split_readsname.txt' % (
        sample, sample))


# 确定split-reads拆分的segment比对到哪些位置
def segment(sample):
    split_read1 = pd.read_csv(os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_split_readsname.txt'),
                              sep='\t', header=None, names=['read', 'cigar'])
    reads_ins = pd.read_csv(os.path.join(workdir, sample, sample + '_ins_reads.txt'),
                            sep='\t', header=None, names=['ins', 'start', 'end', 'read', 'MAPQ', 'strand'])
    new_reads_ins = reads_ins.copy()
    delREAD = []
    split_bed = open(os.path.join(workdir, sample, sample + '_ins_split-reads_detail.bed'), 'w')
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
    new_reads_ins.to_csv(os.path.join(workdir, sample, sample + '_ins_reads_nonsplit.txt'), sep='\t', header=None,
                         index=None)
    os.system(
        'sort -u %s | sort -k1,1 -k2,2n > %s' % (os.path.join(workdir, sample, sample + '_ins_split-reads_detail.bed'),
                                                 os.path.join(workdir, sample,
                                                              'sorted_uniq_' + sample + '_ins_split-reads_detail.bed')))


def new_mapping(sample):
    subworkdir = sample_list[sample]
    new_subworkdir = os.path.join(workdir, sample)
    if not os.path.exists(new_subworkdir):
        os.mkdir(new_subworkdir)
    os.chdir(new_subworkdir)
    fq1 = os.path.join(subworkdir, 'merge_%s_unmap_1.fq' % sample)
    fq2 = os.path.join(subworkdir, 'merge_%s_unmap_2.fq' % sample)
    newfasta = sample + '_all_extendINS.fasta'
    new_fasta = open(newfasta, 'w')
    if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
        Sample = sample
    else:
        Sample = sample + 'C'
    FASTAdir = {'INS': os.path.join(ins_locdir, Sample),
                'DUP': os.path.join(dup_locdir, Sample)}
    for sv in FASTAdir:
        fastadir = FASTAdir[sv]
        if os.path.exists(fastadir):
            for refname in os.listdir(fastadir):
                ins_extend = os.path.join(ins_locdir, Sample, refname, refname + '_SRM.fasta')
                locfile = os.path.join(ins_locdir, Sample, refname, refname + '_SRM_INS_coordinate.bed')
                if os.path.exists(ins_extend) and os.path.exists(locfile):
                    if not any(os.popen("grep '\-1' %s" % locfile).read()):
                        os.system('awk -F \'\t\' \'{print "%s\t"$2"\t"$3}\' %s >> %s_ins_v2n.bed' % (
                        refname, locfile, sample))
                        n = 1
                        for line in open(ins_extend).readlines():
                            if ">" in line:
                                if n == 1:
                                    newline = '>' + refname + '\n'
                                    n += 1
                                elif n > 1:
                                    newline = '>' + refname + '-' + str(n) + '\n'
                                    n += 1
                                new_fasta.write(newline)
                            else:
                                new_fasta.write(line)
    new_fasta.close()
    os.system(hisat2_build + " %s %s -p 10" % (newfasta, new_subworkdir))
    sam = sample + '_unmapReads_mapped_20210624.sam'
    os.system(hisat2 + ' -p 10 -x {sample_workdir} -1 {fq1} -2 {fq2} -S {sam}'.format(
        sample_workdir=new_subworkdir, fq1=fq1, fq2=fq2, sam=sam))
    os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (sam, sam.split(".sam")[0]))
    os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (sam.split(".sam")[0], sam.split(".sam")[0]))
    # os.system('rm *sam')
    split_reads(sample, 'sorted_%s.bam' % sam.split(".sam")[0])  # 提取bam中split reads
    os.system(bedtools + ' bamtobed -i sorted_{idx}.bam > sorted_{idx}.bed'.format(idx=sam.split(".sam")[0]))
    os.system(
        bedtools + ' intersect -a sorted_{idx}.bed -b {sample}_ins_v2n.bed -wa | sort -u > {sample}_ins_reads.txt'.format(
            idx=sam.split(".sam")[0], sample=sample))
    print(sample + ' start')
    segment(sample)  # ins中split reads比对到组装INS什么位置
    os.system(bedtools + ' intersect -a sorted_uniq_{sample}_ins_split-reads_detail.bed \
                    -b {sample}_ins_v2n.bed -wa | sort -u > \
                    {sample}_ins_split-reads_detail.txt'.format(sample=sample))  # 有segment落在INS内的split-read
    print(sample + ' end')


pools = Pool(5)
for sample in list(sample_list.keys()):
    pools.apply_async(new_mapping, args=(sample,))

pools.close()
pools.join()
del pools

# 统计INS内RNA-seq reads数量
reads_result = open('samples_reads_INS_stat_20210702.txt', 'w')
reads_result.write('sample\tchr\tstart\tend\tid\treads\tsplit_reads\n')
for sample in list(sample_list):
    if ('Lung' in sample) and ("ht2" not in sample):
        # os.system(samtools + ' flagstat {workdir}/{sample}/sorted_{sample}_unmapReads_mapped_20210421.bam > \
        # {workdir}/{sample}/sorted_{sample}_unmapReads_mapped_20210421.flagstat'.format(
        #    sample=sample,workdir=workdir))
        reads = os.path.join(workdir, sample,
                             sample + '_ins_reads_nonsplit.txt')
        v2nBed = os.path.join(workdir, sample,
                              sample + '_ins_v2n.bed')
        v2nBed_break = os.path.join(workdir, sample,
                                    sample + '_ins_v2n_break.bed')
        os.system("awk -F '\t' '{print $1\"\t\"$2-10\"\t\"$2+10}' %s >> %s" % (v2nBed, v2nBed_break))
        os.system("awk -F '\t' '{print $1\"\t\"$3-10\"\t\"$3+10}' %s >> %s" % (v2nBed, v2nBed_break))
        newReads = os.path.join(workdir, sample,
                                sample + '_ins_reads_nonsplit_in_INS.txt')
        split = os.path.join(workdir, sample,
                             sample + '_ins_split-reads_detail.txt')
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

out = open('insExtend_site_samples_reads_INS_stat_20210702.bed', 'w')  # 要求nonsplit read必须在INS之内或者跨过
for line in open("samples_reads_INS_stat_20210702.txt").readlines():
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
        hg38_loc = open(os.path.join(ins_locdir, Sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')).readline()
        if '-1' not in hg38_loc:
            newline = hg38_loc.strip() + '\t' + '\t'.join(info[0:7]) + '\n'
            out.write(newline)

out.close()

genebed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_gene_20200312.bed'
gene5kb_bed = '/data/fs08/wangzf/nanopore/xl/ref/gencode_v33_annotation_gtf_gene5kb_cosmic-ncg6.bed'

os.system('sort -k1,1 -k2,2n insExtend_site_samples_reads_INS_stat_20210702.bed > \
            sorted_samples_reads_INS_stat_2021702.bed')
os.system(bedtools + ' closest -a sorted_samples_reads_INS_stat_2021702.bed \
            -b %s -t all -D ref > gene_samples_reads_INS_stat_20210702.bed' % genebed)
os.system(bedtools + ' intersect -a sorted_samples_reads_INS_stat_2021702.bed \
            -b %s -wa -wb | sort -k1,1 -k2,2n -u > gene5kb_overlap_samples_reads_INS_stat_20210702.bed' % gene5kb_bed)

## 确定split segment 比对位置在ref的位置
genebed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_gene_20200312.bed'
exonbed = '/data/fs08/wangzf/hg38_ztf/gencode.v33.annotation_exon_20200312.bed'
bedtools = '/data/fs01/biosoft/bedtools-2.28.0/bin/bedtools'


def segment_gene(sample):
    print(sample)
    if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
        Sample = sample
    else:
        Sample = sample + 'C'
    INSsplit_bed = pd.read_csv(os.path.join(workdir, sample,
                                            sample + '_ins_split-reads_detail.txt'),
                               sep='\t', header=None, names=['refname', 'start', 'end', 'read'])
    ALLsplit_bed = pd.read_csv(os.path.join(workdir, sample, 'sorted_uniq_' + sample + '_ins_split-reads_detail.bed'),
                               sep='\t', header=None, names=['refname', 'start', 'end', 'read'])
    split_bed = ALLsplit_bed[ALLsplit_bed['read'].isin(list(INSsplit_bed['read']))]
    new_split_bed = open(os.path.join(workdir, sample,
                                      'hg38_' + sample + '_ins_split-reads_detail_20210507.bed'), 'w')
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
                                      '_ins_split-reads_detail_20210507.bed')
    New_split_bed = os.path.join(workdir, sample, 'hg38_' + sample +
                                 '_ins_split-reads_detail_20210507.bed')
    exon_new_split_bed = os.path.join(workdir, sample, "exon_hg38_" + sample +
                                      '_ins_split-reads_detail_20210507.bed')
    os.system(bedtools + ' intersect -a %s -b %s -wa -wb > %s' % (New_split_bed, genebed, gene_new_split_bed))
    os.system(bedtools + ' intersect -a %s -b %s -wa -wb > %s' % (New_split_bed, exonbed, exon_new_split_bed))


pools = Pool(5)
for sample in list(sample_list.keys()):
    pools.apply_async(segment_gene, args=(sample,))

pools.close()
pools.join()
del pools


## 确定比对上的基因与INS之间的距离
def ins_gene_exon_distance(sample):
    gene_new_split_bed = os.path.join(workdir, sample, "gene_hg38_" + sample +
                                      '_ins_split-reads_detail_20210507.bed')
    exon_new_split_bed = os.path.join(workdir, sample, "exon_hg38_" + sample +
                                      '_ins_split-reads_detail_20210507.bed')
    file_list = {'gene': gene_new_split_bed, 'exon': exon_new_split_bed}
    if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
        Sample = sample
    else:
        Sample = sample + 'C'
    for fileIDX in file_list:
        file = open(file_list[fileIDX])
        newfile = open(file_list[fileIDX].split(".bed")[0] + '_distance.tmp', 'w')
        for line in file.readlines():
            info = line.strip().split('\t')
            refname = info[3]
            idxstart = int(info[8])
            idxend = int(info[9])
            loc_hg38 = os.path.join(ins_locdir, Sample, refname,
                                    refname + '_SRM_INS_hg38_coordinate.bed')
            if not os.path.exists(loc_hg38):
                loc_hg38 = os.path.join(dup_locdir, Sample, refname,
                                        refname + '_SRM_INS_hg38_coordinate.bed')
            hg38_loc = pd.read_csv(loc_hg38, sep='\t', header=None)
            hg38Start = hg38_loc.iloc[0, 1]
            hg38End = hg38_loc.iloc[0, 2]
            if (idxend < hg38Start):
                newline = line.strip() + '\t' + str(hg38Start - idxend) + '\n'
            elif (idxstart > hg38End):
                newline = line.strip() + '\t' + str(idxstart - hg38End) + '\n'
            else:
                newline = line.strip() + '\t0\n'
            newfile.write(newline)
        newfile.close()
        os.system('sort -u -k1,1 -k2,2n {idx}_distance.tmp > {idx}_distance.bed'.format(
            idx=file_list[fileIDX].split(".bed")[0]))
        os.system('rm {idx}_distance.tmp'.format(idx=file_list[fileIDX].split(".bed")[0]))


pools = Pool(5)
for sample in list(sample_list.keys()):
    pools.apply_async(ins_gene_exon_distance, args=(sample,))

pools.close()
pools.join()
del pools

## 统计INS的类别
dup_type = pd.read_csv("DUP_repeatType_20210624.bed", header=None, sep='\t')
dup_type_all = pd.read_csv("DUP_repeat_20210624.bed", header=None, sep='\t')
ins_type = pd.read_csv("INS_repeatType_20210624.bed", header=None, sep='\t')
ins_type_all = pd.read_csv("INS_repeat_20210624.bed", header=None, sep='\t')


def ins_te(sample):
    gene_new_split_bed = os.path.join(workdir, sample, "gene_hg38_" + sample +
                                      '_ins_split-reads_detail_20210507_distance.bed')
    exon_new_split_bed = os.path.join(workdir, sample, "exon_hg38_" + sample +
                                      '_ins_split-reads_detail_20210507_distance.bed')
    file_list = {'gene': gene_new_split_bed, 'exon': exon_new_split_bed}
    for fileIDX in file_list:
        gene = open(file_list[fileIDX])
        te_gene = open(os.path.join(workdir, sample, "repeat_te_" + fileIDX + "_hg38_" + sample +
                                    '_ins_split-reads_detail_20210507_distance.bed'), 'w')
        for line in gene.readlines():
            info = line.strip().split('\t')
            if info[3] in list(dup_type_all[0]):
                info.extend(dup_type_all[dup_type_all[0] == info[3]].iloc[0].tolist()[3:5])
            else:
                info.extend(ins_type_all[ins_type_all[0] == info[3]].iloc[0].tolist()[3:5])
            repeat = info[-2].split("|")
            rep_list = []
            for rep in repeat:
                rep_list.append('_'.join(rep.split("_")[2:len(rep.split("_"))]))
            info[-2] = '|'.join(rep_list)
            if len(set(rep_list)) == 1:
                newline = '\t'.join(info) + '\ttype_single\n'
            else:
                newline = '\t'.join(info) + '\ttype_multiple\n'
            te_gene.write(newline)
        gene.close()
        te_gene.close()


pools = Pool(5)
for sample in list(sample_list.keys()):
    pools.apply_async(ins_te, args=(sample,))

pools.close()
pools.join()
del pools

## 统计共享情况
patient = pd.read_csv("/data/fs08/wangzf/nanopore/xl/ref/LUAD_Patients_ONT20210120_FORwanguo.csv", sep=',')
patient.iloc[87, 1] = 'Lung61C'

os.system('cat */repeat_te_gene_hg38_*_ins_split-reads_detail_20210507_distance.bed > \
                allSamples_repeat_gene_hg38_split-reads_detail.txt')
all_repeat_gene = pd.read_csv("allSamples_repeat_gene_hg38_split-reads_detail.txt", header=None, sep='\t')
all_repeat_gene1 = all_repeat_gene[[3, 6, 12, 13, 14, 15, 16, 17]]
all_repeat_gene1.columns = ['svid', 'readid', 'genetype', 'gene', 'distance', 'repeat', 'rep_class', 'type_class']
all_repeat_gene1['sample'] = all_repeat_gene1['svid'].str.split("_").str[0]
all_repeat_gene2 = all_repeat_gene1.drop_duplicates()
all_repeat_gene3 = all_repeat_gene2.pivot_table(index=['gene', 'sample'],
                                                values='readid',
                                                aggfunc={'readid': 'count'},
                                                fill_value=0)
all_repeat_gene3['gene'] = all_repeat_gene3.index.get_level_values('gene')
all_repeat_gene3['sample'] = all_repeat_gene3.index.get_level_values('sample')
all_repeat_gene3 = all_repeat_gene3[['gene', 'sample', 'readid']]
all_repeat_gene3.index = list(range(0, all_repeat_gene3.shape[0]))
all_repeat_gene_type = all_repeat_gene1[['gene', 'sample', 'distance', 'repeat', 'rep_class', 'type_class']]
all_repeat_gene_type = all_repeat_gene_type.drop_duplicates()
all_repeat_gene4 = pd.merge(all_repeat_gene3, all_repeat_gene_type, on=['gene', 'sample'])
all_repeat_gene4.to_csv("allSamples_repeat_gene_hg38_split-reads_detail_stat.csv", index=0)

colnames = ['gene', 'genetype', 'sample_count', 'samples', 'cancer_type', 'TNM', 'SV', 'TEfamily', 'split-gene']
gene5kb = pd.read_csv("gene5kb_overlap_samples_reads_INS_stat_20210702.bed", sep='\t', header=None)
gene5kb1 = gene5kb[[3, 4, 5, 6, 7, 8, 9, 15, 16]]
gene5kb1.columns = ['sample', 'chr', 'start', 'end', 'id', 'read', 'split-read', 'genetype', 'gene']
gene5kb1['newid_xl'] = gene5kb1['sample'] + '_' + gene5kb1['chr'] + '_' + gene5kb1['start'].map(str) + \
                       '_' + gene5kb1['end'].map(str) + '_' + gene5kb1['id'].map(str)
gene5kb1 = gene5kb1[['newid_xl', 'read', 'split-read', 'genetype', 'gene']]
gene5kb2 = pd.merge(gene5kb1, teinfo1, on=['newid_xl'], how='left')
gene5kb2 = gene5kb2[gene5kb2['decision_type'].notna()]
gene5kb2_1 = gene5kb2[['genetype', 'gene', 'sample_id']]
gene5kb2_1 = gene5kb2_1.drop_duplicates()
gene5kb3 = gene5kb2_1.pivot_table(index=['genetype', 'gene'],
                                  values='sample_id',
                                  aggfunc={'sample_id': 'count'},
                                  fill_value=0)
gene5kb3['gene'] = gene5kb3.index.get_level_values('gene')
gene5kb3['genetype'] = gene5kb3.index.get_level_values('genetype')
gene5kb3.index = list(range(0, gene5kb3.shape[0]))
stat_df = pd.DataFrame({'gene': [], 'genetype': [], 'sample_count': [],
                        'samples': [], 'cancer_type': [], 'TNM': [],
                        'SV': [], 'TEfamily': [], 'class': [], 'split-gene': []})
for i in range(0, gene5kb3.shape[0]):
    info = gene5kb3.iloc[i].to_list()
    subdf = gene5kb2[gene5kb2['gene'] == info[1]]
    info_list = {'gene': [info[1]], 'genetype': [info[2]], 'sample_count': [info[0]],
                 'samples': [';'.join(list(set(subdf['sample_id'])))],
                 'cancer_type': [], 'TNM': [],
                 'SV': [], 'TEfamily': [],
                 'split-gene': [';'.join(list(all_repeat_gene3[all_repeat_gene3['gene'] == info[1]]['sample']))]}
    cancer_type = []
    TNM = []
    SV = []
    TE = []
    Class = []
    for j in range(0, subdf.shape[0]):
        newinfo = subdf.iloc[j].to_list()
        if (newinfo[-1] == 1):
            sv = 'repeat_dup'
        elif (newinfo[-2] == 1):
            sv = 'segment_dup'
        else:
            sv = 'INS'
        sample = newinfo[7]
        cancer_type.append(patient[patient['sample'] == sample].values.tolist()[0][-5])
        TNM.append(patient[patient['sample'] == sample].values.tolist()[0][-1])
        SV.append(sv)
        TE.append(newinfo[-3])
        repeat = newinfo[-3].split(",")
        if len(set(repeat)) == 1:
            Class.append('type_single')
        else:
            Class.append('type_multiple')
    info_list['TEfamily'] = ';'.join(TE)
    info_list['TNM'] = ';'.join(TNM)
    info_list['cancer_type'] = ';'.join(cancer_type)
    info_list['SV'] = ';'.join(SV)
    info_list['class'] = ';'.join(Class)
    info_list = pd.DataFrame(info_list)
    stat_df = stat_df.append(info_list)

stat_df.to_csv("stat_gene5kb_overlap_samples_reads_INS_stat_20210702.csv", index=0)

################################ stat
os.chdir('/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524/INS_DUP_RNA_20210624')
gene5kb = pd.read_csv("gene5kb_overlap_samples_reads_INS_stat_20210702.bed", sep='\t', header=None)
gene5kb1 = gene5kb[[3, 4, 5, 6, 7, 8, 9, 15, 16]]
gene5kb1.columns = ['sample', 'chr', 'start', 'end', 'id', 'read', 'split-read', 'genetype', 'gene']
gene5kb1['newid_xl'] = gene5kb1['sample'] + '_' + gene5kb1['chr'] + '_' + gene5kb1['start'].map(str) + \
                       '_' + gene5kb1['end'].map(str) + '_' + gene5kb1['id'].map(str)
gene5kb2 = gene5kb1[['gene', 'genetype', 'newid_xl', 'read', 'split-read']]

all_repeat_gene = pd.read_csv("allSamples_repeat_gene_hg38_split-reads_detail.txt", header=None, sep='\t')
all_repeat_gene1 = all_repeat_gene[[3, 6, 12, 13, 14, 15, 16, 17]]
all_repeat_gene1.columns = ['svid', 'readid', 'genetype', 'gene', 'distance', 'repeat', 'rep_class', 'type_class']
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
all_repeat_gene3.columns = ['gene', 'genetype', 'newid_xl', 'split_read_gene']
all_repeat_gene3.index = list(range(0, all_repeat_gene3.shape[0]))

newdf = pd.merge(gene5kb2, all_repeat_gene3, on=['gene', 'genetype', 'newid_xl'], how='outer')
newdf1 = newdf.fillna(0)
newdf1['sample'] = newdf1['newid_xl'].str.split("_").str[0]

dup_type_all = pd.read_csv("DUP_repeat_20210624.bed", header=None, sep='\t')
ins_type_all = pd.read_csv("INS_repeat_20210624.bed", header=None, sep='\t')
outfile = open('allSamples_repeat_gene_INS_stat_20210702.txt', 'w')
outfile.write('#' + '\t'.join(['chr', 'start', 'end', 'gene', 'genetype',
                               'newid_xl', 'read', 'split-read', 'split_read_gene', 'sample',
                               'repeat', 'repeat_class', 'repeatType']) + '\n')
for i in range(0, newdf1.shape[0]):
    info = newdf1.iloc[i].map(str).to_list()
    refname = info[2]
    sample = info[-1]
    loc_hg38 = os.path.join(ins_locdir, sample, refname,
                            refname + '_SRM_INS_hg38_coordinate.bed')
    if not os.path.exists(loc_hg38):
        loc_hg38 = os.path.join(dup_locdir, sample, refname,
                                refname + '_SRM_INS_hg38_coordinate.bed')
    if refname in list(ins_type_all[0]):
        te = ins_type_all[ins_type_all[0] == refname].iloc[0].map(str).to_list()[3:5]
    elif refname in list(dup_type_all[0]):
        te = dup_type_all[dup_type_all[0] == refname].iloc[0].map(str).to_list()[3:5]
    else:
        te = ['unmasked', 'single']
    hg38loc = pd.read_csv(loc_hg38, sep='\t', header=None)
    hg38_loc = hg38loc.iloc[0].map(str).to_list()
    hg38_loc.extend(info)
    info = hg38_loc
    info.extend(te)
    te_list = te[0].split("|")
    new_te_list = []
    for TE in te_list:
        if 'unmasked' not in TE:
            TE = '_'.join(TE.split("_")[2:len(TE.split("_"))])
        new_te_list.append(TE)
    if len(set(new_te_list)) > 1:
        info.append('type_multiple')
    else:
        info.append(new_te_list[0])
    newline = '\t'.join(info) + '\n'
    outfile.write(newline)

outfile.close()

anno = {'gene': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_gene_20200312.bed',
        'exon': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_exon_20200312.bed',
        'promoter': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_promoter_20200312.bed',
        'intron': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_intron_20200312.bed',
        'downstream': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_gene_downstream2kb_20210318.bed',
        'intergenic': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_intergenic_20200312.bed',
        'enhancer': '/data/fs08/wangzf/nanopore/xl/ref/human_permissive_enhancers_phase_1_and_2_LIFTOVER_hg38_new.bed',
        'openPeak': '/data/fs08/wangzf/Assembly_ztf/ATAC_peak/inter_samples_summits_merged_filtered.bed'}

allINS_rna = 'allSamples_repeat_gene_INS_stat_20210702.txt'
for dis in [100, 200, 500]:
    all = pd.read_csv(allINS_rna, sep='\t')
    newbed = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/allSamples_repeat_gene_INS_stat_' + \
             str(dis) + '_20210702.bed'
    all['start'] = all['start'] - dis
    all['end'] = all['end'] + dis
    all.to_csv(newbed, sep='\t', index=0, header=0)

total_ins_rna = '1576'
outfile = open('allSamples_repeat_gene_INS_stat_for-Oddratio_20210702.txt', 'w')
for annoregion in anno:
    os.system(bedtools + ' intersect -a %s -b %s -wa -wb > \
                       /data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/%s_%s' %
              (allINS_rna, anno[annoregion], annoregion, allINS_rna))
    count = os.popen("awk -F '\t' '{print $6}' /data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/%s_%s \
                        | sort -u | wc -l " % (annoregion, allINS_rna)).read().strip().split()[0]
    newline = '\t'.join([annoregion, '0', total_ins_rna, count]) + '\n'
    outfile.write(newline)
    for dis in [100, 200, 500]:
        anno_result = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/allSamples_repeat_gene_INS_stat_' + \
                      str(dis) + '_20210702.bed'
        os.system(bedtools + ' intersect -a %s -b %s -wa -wb > \
                                /data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/%s_%s' % (
            anno_result, anno[annoregion], annoregion, anno_result.split("/")[-1]))
        count = os.popen("awk -F '\t' '{print $6}' /data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/unmapped_reads/%s_%s \
                                | sort -u | wc -l " % (annoregion, anno_result.split("/")[-1])).read().strip().split()[
            0]
        newline = '\t'.join([annoregion, str(dis), total_ins_rna, count]) + '\n'
        outfile.write(newline)

outfile.close()

######### all RNA reads
hg38_te = '/data/fs08/wangzf/nanopore/pxn/02.ref/hg38.repeatmasker.reannot.sorted.bed'
rna_fq_dir = '/data/fs02/wangzf/nanopore_wxy/RNA-seq/cleandata'
workdir = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630'
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
    new_subworkdir = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524/INS_DUP_RNA_20210624/%s' % sample
    sam = sample + '_discReads_mapped_20210630.sam'
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
    bam = '/data/fs08/wangzf/nanopore/rna/RNA-pipeline_bam/%s.sorted.bam' % sample
    pools.apply_async(preparation, args=(sample, bam,))

pools.close()
pools.join()
del pools

seqtk = '/data/fs08/wangzf/nanopore/xl/tools/seqtk/seqtk'
workdir = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630'
exonBed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_exon_20200312.bed'
rna_fq_dir = '/data/fs02/wangzf/nanopore_wxy/RNA-seq/cleandata'


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
    new_subworkdir = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524/INS_DUP_RNA_20210624/%s' % sample
    sam = sample + '_splitReads_mapped_20210630.sam'
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
genebed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_gene_20200312.bed'
gene5kb_bed = '/data/fs08/wangzf/nanopore/xl/ref/gencode_v33_annotation_gtf_gene5kb_cosmic-ncg6.bed'
oldworkdir = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524//INS_DUP_RNA_20210624'

for idx in ['all-split', 'disc']:
    reads_result = open('samples_reads_INS_stat_%s_20210704.txt' % idx, 'w')
    reads_result.write('sample\tchr\tstart\tend\tid\treads\tsplit_reads\n')
    for sample in list(sample_list):
        if ('Lung' in sample) and ("ht2" not in sample):
            # os.system(samtools + ' flagstat {workdir}/{sample}/sorted_{sample}_unmapReads_mapped_20210421.bam > \
            # {workdir}/{sample}/sorted_{sample}_unmapReads_mapped_20210421.flagstat'.format(
            #    sample=sample,workdir=workdir))
            reads = os.path.join(workdir, sample,
                                 sample + '_ins_reads.%s_nonsplit.txt' % idx)
            v2nBed_break = os.path.join(oldworkdir, sample,
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
    out = open('insExtend_site_samples_reads.%s_INS_stat_20210704.bed' % idx, 'w')  # 要求nonsplit read必须在INS之内或者跨过
    for line in open("samples_reads_INS_stat_%s_20210704.txt" % idx).readlines():
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
    os.system('sort -k1,1 -k2,2n insExtend_site_samples_reads.%s_INS_stat_20210704.bed > \
                sorted_samples_reads.%s_INS_stat_2021704.bed' % (idx, idx))
    os.system(bedtools + ' closest -a sorted_samples_reads.%s_INS_stat_2021704.bed \
                -b %s -t all -D ref > gene_samples_reads.%s_INS_stat_20210704.bed' % (idx, genebed, idx))
    os.system(bedtools + ' intersect -a sorted_samples_reads.%s_INS_stat_2021704.bed \
                -b %s -wa -wb | sort -k1,1 -k2,2n -u > \
                gene5kb_overlap_samples_reads.%s_INS_stat_20210704.bed' % (idx, gene5kb_bed, idx))

## 确定split segment 比对位置在ref的位置
genebed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_gene_20200312.bed'
exonbed = '/data/fs08/wangzf/hg38_ztf/gencode.v33.annotation_exon_20200312.bed'
bedtools = '/data/fs01/biosoft/bedtools-2.28.0/bin/bedtools'


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
                                          'hg38_' + sample + '_ins_split-reads.%s_detail_20210704.bed' % idx), 'w')
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
                                          '_ins_split-reads.%s_detail_20210704.bed' % idx)
        New_split_bed = os.path.join(workdir, sample, 'hg38_' + sample +
                                     '_ins_split-reads.%s_detail_20210704.bed' % idx)
        exon_new_split_bed = os.path.join(workdir, sample, "exon_hg38_" + sample +
                                          '_ins_split-reads.%s_detail_20210704.bed' % idx)
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
workdir = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630'
ins_locdir = '/data/fs08/wangzf/Assembly_ztf/INS_RepeatMasker_f50bp/Assembly_result'
dup_locdir = '/data/fs08/wangzf/Assembly_DUP/DUP_2286/Assembly_result'
os.chdir(workdir)
oldresult = pd.read_csv(os.path.join(oldworkdir, 'allSamples_repeat_gene_INS_stat_20210702.txt'), sep='\t')
newresult = oldresult.copy()
for idx in ['disc', 'all-split']:
    reads = pd.read_csv('samples_reads_INS_stat_%s_20210704.txt' % idx, sep='\t')
    reads['newid_xl'] = reads['sample'] + '_' + reads['chr'] + '_' + reads['start'].map(str) \
                        + '_' + reads['end'].map(str) + '_' + reads['id'].map(str)
    reads.columns = ['sample', 'chr', 'start', 'end', 'id', 'disc-reads', 'disc-split_reads', 'newid_xl']
    gene5kb = pd.read_csv("gene5kb_overlap_samples_reads.%s_INS_stat_20210704.bed" % idx, sep='\t', header=None)
    gene5kb1 = gene5kb[[3, 4, 5, 6, 7, 8, 9, 15, 16]]
    gene5kb1.columns = ['sample', 'chr', 'start', 'end', 'id', idx + '_read', idx + '_split-read', 'genetype', 'gene']
    gene5kb1['newid_xl'] = gene5kb1['sample'] + '_' + gene5kb1['chr'] + '_' + gene5kb1['start'].map(str) + \
                           '_' + gene5kb1['end'].map(str) + '_' + gene5kb1['id'].map(str)
    gene5kb2 = gene5kb1[['gene', 'genetype', 'newid_xl', idx + '_read', idx + '_split-read']]
    os.system('cat */gene_hg38_Lung*_ins_split-reads.%s_detail_20210704.bed > \
                    allSamples_gene_hg38_split-reads.%s_detail_20210705.txt' % (idx, idx))
    all_repeat_gene = pd.read_csv("allSamples_gene_hg38_split-reads.%s_detail_20210705.txt" % idx, header=None,
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

dup_type_all = pd.read_csv(oldworkdir + "/DUP_repeat_20210624.bed", header=None, sep='\t')
ins_type_all = pd.read_csv(oldworkdir + "/INS_repeat_20210624.bed", header=None, sep='\t')

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
newresult1.to_csv("allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", sep='\t', index=0)
teinfo2 = teinfo1[['newid_xl', 'segment_dup', 'repeat_dup', 'total_len']]
newresult2 = pd.merge(newresult1, teinfo2, on=['newid_xl'], how='left')
ins_newresult = newresult2[(newresult2['segment_dup'] == 0) & (newresult2['repeat_dup'] == 0)]
dup_newresult = newresult2[(newresult2['segment_dup'] != 0) | (newresult2['repeat_dup'] != 0)]
ins_newresult.to_csv("INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", sep='\t', index=0)
dup_newresult.to_csv("TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", sep='\t', index=0)

anno = {'gene': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_gene_20200312.bed',
        'exon': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_exon_20200312.bed',
        'promoter': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_promoter_20200312.bed',
        'intron': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_intron_20200312.bed',
        'downstream': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_gene_downstream2kb_20210318.bed',
        'intergenic': '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_intergenic_20200312.bed',
        'enhancer': '/data/fs08/wangzf/nanopore/xl/ref/human_permissive_enhancers_phase_1_and_2_LIFTOVER_hg38_new.bed',
        'openPeak': '/data/fs08/wangzf/Assembly_ztf/ATAC_peak/inter_samples_summits_merged_filtered.bed'}

file_list = ['allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt',
             "INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt",
             "TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt"]
for file in file_list:
    for dis in [100, 200, 500]:
        all = pd.read_csv(allINS_rna, sep='\t')
        newbed = 'result_20210705/%s' % file.split('_20210705')[0] + str(dis) + '_20210705.bed'
        all['start'] = all['start'] - dis
        all['end'] = all['end'] + dis
        all.to_csv(newbed, sep='\t', index=0, header=0)

total_ins_rna = '1775'
totalINS = '1423'
totalTD = '352'
new_file_list = {'all': ['allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt', total_ins_rna],
                 'INS': ["INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", totalINS],
                 'TD': ["TD_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", totalTD]}
outfile = open('allSamples_repeat_gene_INS_stat_unmapped_disc_split_for-Oddratio_20210705.txt', 'w')
for annoregion in anno:
    for idx in new_file_list:
        os.system(bedtools + ' intersect -a %s -b %s -wa -wb > result_20210705/%s_%s' %
                  (new_file_list[idx][0], anno[annoregion], annoregion, new_file_list[idx][0]))
        count = os.popen("awk -F '\t' '{print $6}' result_20210705/%s_%s \
                            | sort -u | wc -l " % (annoregion, allINS_rna)).read().strip().split()[0]
        newline = '\t'.join([idx, annoregion, '0', new_file_list[idx][1], count]) + '\n'
        outfile.write(newline)
        for dis in [100, 200, 500]:
            anno_result = 'result_20210705/%s' % new_file_list[idx][0].split('_20210705')[0] + \
                          str(dis) + '_20210705.bed'
            os.system(bedtools + ' intersect -a %s -b %s -wa -wb > result_20210705/%s_%s' % (
                anno_result, anno[annoregion], annoregion, anno_result.split("/")[-1]))
            count = os.popen("awk -F '\t' '{print $6}' result_20210705/%s_%s \
                                    | sort -u | wc -l " % (
            annoregion, anno_result.split("/")[-1])).read().strip().split()[0]
            newline = '\t'.join([idx, annoregion, str(dis), new_file_list[idx][1], count]) + '\n'
            outfile.write(newline)

outfile.close()



#######################################################################################################################
##########
# all RNA reads mapped to candidate INS sequence
gene_list = os.popen("grep -v '#' INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt | \
                        awk -F '\t' '{print $4}' | sort -u").read().strip().split('\n')
insbed = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all_samples_newINS_candidateRegion_20210622.bed'
os.system(bedtools + ' intersect -a %s -b %s -wa -wb > gene5kb_all_samples_newINS_candidateRegion_20210622.txt' % (
gene5kb_bed, insbed))
gene5kb_ins = pd.read_csv("gene5kb_all_samples_newINS_candidateRegion_20210622.txt", sep='\t', header=None)
gene5kb_ins1 = gene5kb_ins[gene5kb_ins[6].isin(gene_list)]
gene5kb_ins1['ins_region'] = gene5kb_ins1[11] + '_' + gene5kb_ins1[12].map(str) + '_' + gene5kb_ins1[13].map(str) +\
                             '_' + gene5kb_ins1[15]
teinfo1['sample'] = teinfo1['newid_xl'].str.split("_").str[0]
teinfo1['new_ins_region'] = teinfo1['ins_region'] + '_' + teinfo1['sample']
ins_teinfo1_gene5kb = teinfo1[teinfo1['new_ins_region'].isin(list(gene5kb_ins1['ins_region']))]
rna_fq_dir = '/data/fs02/wangzf/nanopore_wxy/RNA-seq/cleandata'

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

def bam_png(Bam, Bed, reference, sample, refname, idx):
    chrom, start, end = open(Bed,'r').readline().strip().split('\t')[0:3]
    if os.path.getsize(Bam) < (50 * 1024 * 1024):
        dataset_paths = [Bam, Bed]
        doc = genomeview.visualize_data(dataset_paths, chrom, int(start) - 1000, int(end) + 1000, reference)
        genomeview.save(doc, "_".join([sample, refname, idx]) + ".svg")
        fileHandle = open("_".join([sample, refname, idx]) + ".svg")
        svg = fileHandle.read()
        fileHandle.close()
        exportPath = "_".join([sample, refname, idx]) + ".png"
        exportFileHandle = open(exportPath, 'w')
        cairosvg.svg2png(bytestring=svg, write_to=exportPath)
        exportFileHandle.close()
        os.system("rm %s.svg" % "_".join([A549_sv, refname, idx]))
    else:
        os.system('echo "%s is too large" >> plo.log' % Bam)


pools = Pool(5)
for sample in list(set(ins_teinfo1_gene5kb['sample_id'])):
    subdf = ins_teinfo1_gene5kb[ins_teinfo1_gene5kb['sample_id'] == sample]
    workdir = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/assembled_candidate/'
    for i in range(0, subdf.shape[0]):
        refname = subdf.iloc[i, 0]
        bam = os.path.join(workdir, sample,'sorted_{sample}_allRNAreads_mapped_{refname}.bam'.format(
                            sample=sample, refname=refname))
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        locfile = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
            locfile = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        pools.apply_async(bam_png, args=(Bam, locfile, ins_extend,
                                         sample, refname, 'ins_allreads',))

pools.close()
pools.join()
del pools




##########
# gene structural 对应到 assembled INS上 # 10.100.2.13
ins_rna = pd.read_csv("INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt", sep='\t')
# os.system('grep -v "#" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf > /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_noanno.gtf')
gene_bed = pd.read_csv('/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_noanno.gtf', sep='\t', header=None)
gene_bed['gene'] = gene_bed[8].str.split('gene_name "').str[1]
gene_bed['gene'] = gene_bed['gene'].str.split('"').str[0]
gene_bed['gene_id'] = gene_bed[8].str.split('gene_id "').str[1]
gene_bed['gene_id'] = gene_bed['gene_id'].str.split('"').str[0]

ins_rna1 = ins_rna[['gene', 'newid_xl', 'total_len']]
ins_rna1_gene = pd.merge(ins_rna1, gene_bed, on=['gene'], how='left')
ins_locdir = '/data/fs08/wangzf/Assembly_ztf/INS_RepeatMasker_f50bp/Assembly_result'
dup_locdir = '/data/fs08/wangzf/Assembly_DUP/DUP_2286/Assembly_result'
newgtf = open('INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf', 'w')
exon_newgtf = open('exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf', 'w')
if os.path.exists("newgtf.log"):
    os.system('rm newgtf.log')


for i in range(0, ins_rna1_gene.shape[0]):
    info = ins_rna1_gene.iloc[i].map(str).tolist()
    refname = info[1]
    ins_hg38_bam = os.path.join(ins_locdir, refname.split('_')[0], refname, refname + '_SRM2hg38.bam')
    if not os.path.exists(ins_hg38_bam):
        ins_hg38_bam = os.path.join(dup_locdir, refname.split('_')[0], refname, refname + '_SRM2hg38.bam')
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
os.system('sort -k1,1 -k4,4n -k5,5n INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf > \
            sorted_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf')
os.system('sort -k1,1 -k4,4n -k5,5n exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf > \
            sorted_exon_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf')

############
te_bed = os.path.join(oldworkdir, 'all.INS.DUP.anno.region.clean.20210622.te_list.bed')
genebed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_gene_20200312.bed'

gene_new_split_bed_disc = 'gene_hg38_allSample_ins_split-reads.dics_detail_20210704.bed'
os.system('cat '+workdir+"/*/gene_hg38_*_ins_split-reads.disc_detail_20210704.bed | sort -k1,1 -k2,2n > " +
          gene_new_split_bed_disc)
gene_new_split_bed_all_split = 'gene_hg38_allSample_ins_split-reads.all-split_detail_20210704.bed'
os.system('cat '+workdir+"/*/gene_hg38_*_ins_split-reads.all-split_detail_20210704.bed | sort -k1,1 -k2,2n > " +
          gene_new_split_bed_all_split)
gene_new_split_bed_unmapped = 'gene_hg38_allSample_ins_split-reads.unmapped_detail.bed'
os.system('cat '+oldworkdir+'/*/gene_hg38_*_ins_split-reads_detail_20210507_distance.bed| sort -k1,1 -k2,2n > ' +
          gene_new_split_bed_unmapped)
gene_hg38_reads = {'disc':gene_new_split_bed_disc,
                   'all-split':gene_new_split_bed_all_split,
                   'unmapped':gene_new_split_bed_unmapped}

unmapped_newReads = os.path.join(workdir, 'allSample_ins_reads_nonsplit_in_INS.txt')
os.system('cat %s/*/*_ins_reads_nonsplit_in_INS.txt | sort -k1,1 -k2,2n > %s' % (oldworkdir, unmapped_newReads))
unmapped_split = os.path.join(workdir, 'allSample_ins_split-reads_detail.txt')
os.system('cat %s/*/*_ins_split-reads_detail.txt | sort -k1,1 -k2,2n > %s' % (oldworkdir, unmapped_split))
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

#######################################################################################################################
###################### A549
workdir = '/data/fs09/wangzf/nanopore/xl/A549-new/'
if not os.path.exists(workdir):
    os.mkdir(workdir)


os.chdir(workdir)
# A549 是否携带相同INS
#a549_ontDir = '/data/fs01/wangzf/nanopore/cell_line/A549/'
a549_ontDir = '/data/fs08/wangzf/nanopore/ztf/ascp_A549/A549_DRR171/'

# A549 SV
def sniffles_vcftobed(vcf, bed, tra):
    """vcf to bed
    header: chrom start end (chr start end) svtype id length RE RNAMES IMPRECISE/PRECISE STD_quant_start:STD_quant_stop
    SUPTYPE GT:DR:DV:AF
    """
    bedout = open(bed, 'w')
    traout = open(tra, 'w')
    with open(vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                content = line.split()
                chr1 = content[0]
                start = content[1]
                svid = content[2]
                info = content[7].split(';')
                fil = info[0]
                end = info[3].split('=')[1]
                std1 = info[4].split('=')[1]
                std2 = info[5].split('=')[1]
                svtype = info[8].split('=')[1]
                if '/' in svtype:
                    svtype = '-'.join(svtype.split('/'))
                rnames = info[9].split('=')[1]
                suptype = info[10].split('=')[1]
                svlen = info[11].split('=')[1]
                re = info[14].split('=')[1]
                # af = info[16].split('=')[1]
                std = ':'.join([std1, std2])
                # gt = ':'.join(content[-1].split(':') + [af])
                # anotation = [svtype, svid, svlen, suptype, fil, std, re, rnames, gt]
                anotation = [svtype, svid, svlen, re, rnames, fil, std, suptype]  # , gt]
                if svtype == "TRA":
                    chr2 = info[2].split('=')[1]
                    newline1 = '\t'.join([chr1, start, str(int(start) + 1), chr2, end,
                                          str(int(end) + 1)] + anotation) + '\n'
                    # newline2 = '\t'.join([chr2, end, str(int(end) + 1), chr1, start,
                    #                     str(int(start) + 1)] + anotation) + '\n'
                    # newline = newline1 + newline2
                    traout.write(newline1)
                    traout.flush()
                else:
                    newline = '\t'.join([chr1, start, end] + anotation) + '\n'
                    bedout.write(newline)
                    bedout.flush()
    bedout.close()


def sniffles_svtype(bed, tra, mainbed, maintra):
    """select RE>=5
    delete chrY,chrM
    """
    rmchr = ['chrY', 'chrM']
    mainbedout = open(mainbed, 'w')
    ff = open(bed, 'r')
    line = ff.readline()
    while line:
        l = line.strip().split('\t')
        if l[0] not in rmchr:
            mainbedout.write(line)
            line = ff.readline()
        else:
            line = ff.readline()
    traout = open(maintra, 'w')
    f1 = open(tra, 'r')
    line1 = f1.readline()
    while line1:
        l1 = line1.strip().split('\t')
        if l1[0] in rmchr or l1[3] in rmchr:
            line1 = f1.readline()
        else:
            traout.write(line1)
            line1 = f1.readline()


#vcf = os.path.join(a549_ontDir, 'A549.sorted_ngmlr_sniffles.vcf')
oldvcf = open(os.path.join(a549_ontDir, 'A549_DRR171.sorted_ngmlr_sniffles_s1.vcf'))
vcf = os.path.join(a549_ontDir, 'new_A549_DRR171.sorted_ngmlr_sniffles_s1.vcf')
newvcf = open(vcf,'w')
for line in oldvcf.readlines():
    if '#' not in line:
        info = line.split("\t")
        if info[-2].split(";")[-1]!='RE=1':
            newvcf.write(line)
    else:
        newvcf.write(line)


oldvcf.close()
newvcf.close()
bed = os.path.join(a549_ontDir, 'A549_sv_sniffles_raw.bed')
tra = os.path.join(a549_ontDir, 'A549_tra_sniffles_raw.bed')
sniffles_vcftobed(vcf, bed, tra)
mainbed = os.path.join(a549_ontDir, 'MainChr_A549_sv_sniffles.bed')
maintra = os.path.join(a549_ontDir, 'MainChr_A549_tra_sniffles.bed')
sniffles_svtype(bed, tra, mainbed, maintra)

ins_rna_bed = '/data/fs09/wangzf/nanopore/xl/INS_DUP_RNA_20210630/INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt'
os.system("grep -v '#' %s| awk -F '\t' '{$1\"\t\"$2-500\"\t\"$3+500\"\t\"$4\"\t\"$6\"\t\"$8\"\t\"$9}' | \
            sort -k1,1 -k2,2n > INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.500bp.bed" %
            ins_rna_bed)
os.system(bedtools + ' intersect -a %s -b INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.500bp.bed -wa -wb \
                        | grep -v "INV" | grep -v "DEL" | grep -v "TRA" > \
                        a549_snifflesSV_INS-rna_20210705.txt' % mainbed)
# 提取a549-INS
seqtk = '/data/fs08/wangzf/nanopore/xl/tools/seqtk/seqtk'
ngmlr = "/data/fs01/wangzf/software/ngmlr-0.2.7/ngmlr"
hg38 = '/data/fs01/wangzf/nanopore/ref/mainChr_hg38_UCSC.fa'

def bam_png(chrom, start, end, Bam, Bed, reference, A549_sv, refname, idx):
    if os.path.getsize(Bam) < (50 * 1024 * 1024):
        dataset_paths = [Bam, Bed]
        doc = genomeview.visualize_data(dataset_paths, chrom, int(start) - 1000, int(end) + 1000, reference)
        genomeview.save(doc, "_".join([A549_sv, refname, idx]) + ".svg")
        fileHandle = open("_".join([A549_sv, refname, idx]) + ".svg")
        svg = fileHandle.read()
        fileHandle.close()
        exportPath = "_".join([A549_sv, refname, idx]) + ".png"
        exportFileHandle = open(exportPath, 'w')
        cairosvg.svg2png(bytestring=svg, write_to=exportPath)
        exportFileHandle.close()
        os.system("rm %s.svg" % "_".join([A549_sv, refname, idx]))
    else:
        os.system('echo "%s is too large" >> plo.log' % Bam)


a549_INS_rna = pd.read_csv(workdir + "/a549_snifflesSV_INS-rna_20210705.txt", sep='\t', header=None)
#a549_bam = os.path.join(a549_ontDir, 'A549.sorted_ngmlr.bam')
a549_bam = os.path.join(a549_ontDir, 'A549_DRR171.sorted_ngmlr.bam')
os.system(samtools+' index ' + a549_bam)


def a549_mapped_INS(A549_sv, refname, info):
    subworkdir = os.path.join(workdir, A549_sv + '-' + refname)
    if not os.path.exists(subworkdir):
        os.system('echo "%s------start" >> plot.log' % "_".join([A549_sv, refname, time.asctime()]))
        os.mkdir(subworkdir)
        idx = True
    elif (not os.path.exists(os.path.join(subworkdir, '_'.join([A549_sv, refname, '_hg38.png'])))) or \
            (not os.path.exists(os.path.join(subworkdir, '_'.join([A549_sv, refname, '_INS.png'])))):
        os.system('echo "%s------start" >> plot.log' % "_".join([A549_sv, refname, time.asctime()]))
        idx = True
    else:
        idx = False
    if idx:
        os.chdir(subworkdir)
        sample = refname.split("_")[0]
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        locfile = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
            locfile = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        bam = os.path.join(workdir, '_'.join([A549_sv, 'A549_support_reads.bam']))
        region_bam = os.path.join(workdir, '_'.join([A549_sv, 'A549_support_reads_region.bam']))
        chrom, start, end = A549_sv.split("_")
        Bed = A549_sv + '_sniffles.bed'
        with open(Bed, 'w') as b:
            b.write('\t'.join(info[0:4]) + '\n')
        b.close()
        bam_png(chrom, start, end, region_bam, Bed, hg38, A549_sv, refname, 'hg38')
        support_fq = bam.split(".bam")[0] + '.fq'
        os.system(samtools + ' fastq -@ 10 %s > %s' % (bam, support_fq))
        INSsam = '_'.join([A549_sv, 'A549_support_reads_mappedTo_', refname, '.sam'])
        os.system(
            "{ngmlr} -t 20 -x ont -r {ref} -q {fq} -o {sam}".format(ref=ins_extend, fq=support_fq,
                                                                    sam=INSsam, ngmlr=ngmlr))
        os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (INSsam, INSsam.split(".sam")[0]))
        os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (INSsam.split(".sam")[0], INSsam.split(".sam")[0]))
        os.system(samtools + ' index sorted_%s.bam' % INSsam.split(".sam")[0])
        os.system(
            bedtools + ' bamtobed -split -i sorted_{idx}.bam > sorted_{idx}.bed'.format(idx=INSsam.split(".sam")[0]))
        os.system(
            bedtools + ' intersect -a sorted_{idx}.bed -b {locfile} -wa -wb | sort -u \
                        > {A549_sv}_{refname}_ins_all-split.reads.txt'.format(
                idx=INSsam.split(".sam")[0], A549_sv=A549_sv, refname=refname, locfile=locfile))
        if os.path.getsize('{A549_sv}_{refname}_ins_all-split.reads.txt'.format(A549_sv=A549_sv, refname=refname)) > 0:
            chrom, start, end = os.popen('cat %s' % locfile).read().strip().split("\t")[0:3]
            bam_png(chrom, start, end, "sorted_%s.bam" % INSsam.split(".sam")[0],
                    locfile, ins_extend, A549_sv, refname,
                    'INS')
        os.system('echo "%s------end" >> plot.log' % "_".join([A549_sv, refname, time.asctime()]))


for i in range(0, a549_INS_rna.shape[0]):
    info = a549_INS_rna.iloc[i].map(str).tolist()
    A549_sv = '_'.join(info[0:3])
    reads_list = info[7].split(",")
    refname = info[-3]
    bam = os.path.join(workdir, '_'.join([A549_sv, 'A549_support_reads.bam']))
    region_bam = os.path.join(workdir, '_'.join([A549_sv, 'A549_support_reads_region.bam']))
    if not os.path.exists(region_bam):
        os.system('{samtools} view -@ 10 -b -h {bam} {chrom}:{start}-{end} > {out}'.format(
            samtools=samtools, bam=a549_bam, chrom=A549_sv.split("_")[0],
            start=str(int(A549_sv.split("_")[1]) - 3000),
            end=str(int(A549_sv.split("_")[2]) + 3000),
            out=region_bam))
        os.system("{samtools} index {bam}".format(samtools=samtools, bam=region_bam))
        readBam = pysam.AlignmentFile(region_bam, "rb")
        outbam = pysam.AlignmentFile(bam, 'wb', template=readBam)
        allreads = readBam.fetch()
        n = 0
        for read in allreads:
            if read.query_name in reads_list:
                outbam.write(read)
                n += 1
        readBam.close()
        outbam.close()

pools = Pool(3)
for i in range(0, a549_INS_rna.shape[0]):
    info = a549_INS_rna.iloc[i].map(str).tolist()
    A549_sv = '_'.join(info[0:3])
    reads_list = info[7].split(",")
    refname = info[-3]
    pools.apply_async(a549_mapped_INS, args=(A549_sv, refname, info,))

pools.close()
pools.join()
del pools

a549_genome_support_INS = pd.DataFrame({'sv': [], 'ins': []})
for dir in os.listdir(workdir):
    if ('chr' in dir.split("_")[0]) and os.path.isdir(dir):
        hg38png = os.path.join(workdir, dir, '_'.join(dir.split("-")) + '_hg38.png')
        inspng = os.path.join(workdir, dir, '_'.join(dir.split("-")) + '_INS.png')
        if os.path.exists(hg38png) and os.path.exists(inspng):
            newinfo = dir.split("-")
            a549_genome_support_INS = a549_genome_support_INS.append(
                pd.DataFrame({'sv': [newinfo[0]], 'ins': [newinfo[1]]}))

a549_genome_support_INS['chr'], a549_genome_support_INS['start'], a549_genome_support_INS['end'] = \
    a549_genome_support_INS['sv'].str.split("_").str
a549_genome_support_INS1 = a549_genome_support_INS[['chr', 'start', 'end', 'ins']]
a549_genome_support_INS1['start'] = a549_genome_support_INS1['start'].map(int) - 500
a549_genome_support_INS1['end'] = a549_genome_support_INS1['end'].map(int) + 500
a549_genome_support_INS1 = a549_genome_support_INS1.drop_duplicates()
a549_genome_support_INS1.to_csv("a549_similar_somaticINS.bed", sep='\t', header=0, index=0)
os.mkdir(os.path.join(workdir, 'ONT_NGS_RNA'))

newfasta = './ONT_NGS_RNA/a549_similar_somaticINS_extendINS.fasta'
new_fasta = open(newfasta, 'w')
refname_list = []
for refname in list(a549_genome_support_INS1['ins']):
    if refname not in refname_list:
        refname_list.append(refname)
        sample = refname.split("_")[0]
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        locfile = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
            locfile = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_coordinate.bed')
        if ('C1' in sample) or ('C2' in sample) or ('C3' in sample):
            Sample = sample
        else:
            Sample = sample + 'C'
        if os.path.exists(ins_extend) and os.path.exists(locfile):
            if not any(os.popen("grep '\-1' %s" % locfile).read()):
                os.system(
                    'awk -F \'\t\' \'{print "%s\t"$2"\t"$3}\' %s >> ./ONT_NGS_RNA/a549_similar_somaticINS_ins_v2n.bed' % (
                    refname, locfile))
                n = 1
                for line in open(ins_extend).readlines():
                    if ">" in line:
                        if n == 1:
                            newline = '>' + refname + '\n'
                            n += 1
                        elif n > 1:
                            newline = '>' + refname + '-' + str(n) + '\n'
                            n += 1
                        new_fasta.write(newline)
                    else:
                        new_fasta.write(line)

new_fasta.close()

# NANOPORE RNA: DRX218781: MinION sequencing of SAMD00226192(Aberrant transcript isoforms detected by full-length transcriptome sequencing as transcripts of potential neoantigens in non-small cell lung cancer)
# A549 RNA # 10.100.2.2
a549_ONTrna_fq = '/data/fs09/wangzf/nanopore/xl/A549/DRR228517.fastq'
a549_NGSrna_fq1 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate02_1.fastq.gz'
a549_NGSrna_fq2 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate02_2.fastq.gz'

minimap2 = '/data/fs01/wangzf/software/minimap2/minimap2'
gffcompare = '/data/fs01/wangzf/software/gffcompare-0.12.2.Linux_x86_64/gffcompare'

#os.mkdir(os.path.join(workdir, 'ONT_NGS_RNA'))
os.chdir(os.path.join(workdir, 'ONT_NGS_RNA'))
## ONT hg38
os.system(minimap2 + ' -ax splice -uf -k14 {} {} > A549_ONT_rna.sam'.format(hg38, a549_ONTrna_fq))
os.system(samtools + " view -bS A549_ONT_rna.sam -@ 10 > A549_ONT_rna.bam")
os.system(samtools + " sort -@ 10 A549_ONT_rna.bam > sorted_A549_ONT_rna.bam")
os.system(stringtie_2_1_6 + ' sorted_A549_ONT_rna.bam --conservative -p 10 -o A549_ONT_rna.gtf')
gencode = '/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf'
os.system(gffcompare + ' -R -r %s -o A549_stringtieAssembled_ONT_comp_gencodev33 A549_ONT_rna.gtf' % gencode)
# 提取候选基因相关比较结果
candidate_gene = list(set(a549_INS_rna[14]))  # 525 genes
candidate_gene_ID = list(set(gene_bed[gene_bed['gene'].isin(candidate_gene)]['gene_id']))
tmap = pd.read_csv("A549_stringtieAssembled_ONT_comp_gencodev33.A549_ONT_rna.gtf.tmap", sep='\t')
new_tmap = tmap[tmap['ref_gene_id'].isin(candidate_gene_ID)]
new_tmap.to_csv("filter1_A549_ONT_gencodev33_compare.tmap", sep='\t', index=0)
# 过滤，只保留exon>1并且长度>200bp的transcripts
os.system("awk '($6>1 && $10>=200){print$0}' filter1_A549_ONT_gencodev33_compare.tmap \
            > filter2_by_exon_length_A549_ONT.tmap")
os.system("awk '{print $5}' filter2_by_exon_length_A549_ONT.tmap > filter2_by_exon_length_A549_ONT_transcript_ID")
os.system("grep -Ff filter2_by_exon_length_A549_ONT_transcript_ID -w A549_ONT_rna.gtf > \
            filter2_by_exon_length_A549_ONT_transcript.gtf")
filter2_tmap = pd.read_csv("filter2_by_exon_length_A549_ONT.tmap", sep='\t')

## ONT-a549 INS
os.system(minimap2 + ' -ax splice -uf -k14 a549_similar_somaticINS_extendINS.fasta {} >\
                        A549_ONT_rna-a549_similar_INS.sam'.format(a549_ONTrna_fq))
os.system(samtools + " view -bS A549_ONT_rna-a549_similar_INS.sam -@ 10 > A549_ONT_rna-a549_similar_INS.bam")
os.system(samtools + " sort A549_ONT_rna-a549_similar_INS.bam > sorted_A549_ONT_rna-a549_similar_INS.bam")
os.system(samtools + ' index sorted_A549_ONT_rna-a549_similar_INS.bam')
os.system(bedtools + ' bamtobed -split -i sorted_A549_ONT_rna-a549_similar_INS.bam | sort -k1,1 -k2,2n > \
                        sorted_A549_ONT_rna-a549_similar_INS.bed')
os.system('sort -k1,1 -k2,2n a549_similar_somaticINS_ins_v2n.bed > sorted_a549_similar_somaticINS_ins_v2n.bed')
# if ONT reads mapped on INS seq overlapped INS region
os.system(bedtools + ' intersect -a sorted_A549_ONT_rna-a549_similar_INS.bed \
                        -b a549_similar_somaticINS_ins_v2n.bed -wa -wb > \
                        targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bed') # read bed overlapped with INS regions
target_INS = pd.read_csv("targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bed",sep='\t',header=None,
                         usecols = [0,1,2,3],names=['refname','start','end','read'])

os.system('/data/fs01/biosoft/sambamba-0.7.0/sambamba slice -L a549_similar_somaticINS_ins_v2n.bed \
            sorted_A549_ONT_rna-a549_similar_INS.bam > targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bam')
reads = list(set(target_INS['read']))
samfile = pysam.AlignmentFile("sorted_A549_ONT_rna-a549_similar_INS.bam", "rb")
outbam =  pysam.AlignmentFile("a549_ont_support_RNA_reads_stat_genename-readsName.bam", "wb", template=samfile)
allreads=samfile.fetch()
n=0
for read in allreads:
    if read.query_name in reads:
        outbam.write(read)
        n +=1


samfile.close()
outbam.close()
os.system(samtools + ' index a549_ont_support_RNA_reads_stat_genename-readsName.bam')

os.system(stringtie + " sorted_A549_ONT_rna-a549_similar_INS.bam -m 100 -c 1 -o sorted_A549_ONT_rna-a549_similar_INS.gtf -p 10")
os.system(gffcompare + ' -R -r %s -o  sorted_A549_ONT_rna-a549_similar_INS_comp_gencodev33 sorted_A549_ONT_rna-a549_similar_INS.gtf' % gencode)


# if ONT reads mapped on INS region mapped on known exon
newgtf = '../../INS_DUP_RNA_20210630/INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf'
newgtf_df = pd.read_csv(newgtf, sep='\t',header=None)
similar_gtf = 'new_similar_INS-rna_gencodeV33.gtf'
ins_rna_df = pd.read_csv("../a549_similar_somaticINS.bed", sep='\t',header=None)
newgtf_df1 = newgtf_df[newgtf_df[0].isin(list(ins_rna_df[3]))]
newgtf_df1.to_csv(similar_gtf,index=0,header=0,sep='\t')

gtf = 'new_similar_INS-rna_gencodeV33.gtf'
gtf_df = pd.read_csv(gtf, sep='\t',header=None)
exon_gtf = gtf_df[gtf_df[2]=='exon']
exon_gtf = exon_gtf[[0,3,4,2,8]]
exon_gtf['gene'] = exon_gtf[8].str.split('gene_name "').str[1]
exon_gtf['gene'] = exon_gtf['gene'].str.split('"').str[0]
exon_gtf['genetype'] = exon_gtf[8].str.split('gene_type "').str[1]
exon_gtf['genetype'] = exon_gtf['genetype'].str.split('"').str[0]
exon_gtf['transcript_id'] = exon_gtf[8].str.split('transcript_id "').str[1]
exon_gtf['transcript_id'] = exon_gtf['transcript_id'].str.split('"').str[0]
exon_gtf['exon_number'] = exon_gtf[8].str.split('exon_number ').str[1]
exon_gtf['exon_number'] = exon_gtf['exon_number'].str.split(';').str[0]
exon_gtf=exon_gtf[[0,3,4,2,'gene','genetype','transcript_id','exon_number']]
exon_gtf.to_csv("new_similar_INS-rna_gencodeV33.exon.bed",sep='\t',index=0,header=0)

ins_reads = pd.read_csv("sorted_A549_ONT_rna-a549_similar_INS.bed",sep='\t',header=None)
ins_region_reads = pd.read_csv("targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bed",sep='\t',header=None)
new_ins_reads = ins_reads[ins_reads[3].isin(list(set(ins_region_reads[3])))]
new_ins_reads.to_csv('targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed',sep='\t',header=0,index=0)
os.system(bedtools +' intersect -a targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed \
                        -b new_similar_INS-rna_gencodeV33.exon.bed -wa -wb >\
                         exon_targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed')
exon_ins_region_reads = pd.read_csv("exon_targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed",sep='\t',header=None)
ins_ont_stat = pd.DataFrame({"chr":[],"start":[],"end":[],"INS":[],'INS_type':[],
                             'support_read_count':[],'gene':[],'genetype':[],
                             'exon_count':[],'exon_count_detail':[],
                             'sample_shared':[],'sample_gene':[],'sample_genetype':[]})
ins_bed = pd.read_csv('/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all.INS.DUP.anno.region.clean.20210524.tsv',
                      sep='\t',header=0)
oldworkdir = '/data/fs08/wangzf/nanopore/xl/new_somatic_SV_TGS_bed_20210524//INS_DUP_RNA_20210624'
ins_type_all = pd.read_csv(oldworkdir + "/INS_repeat_20210624.bed", header=None, sep='\t')
ins_genestat = pd.read_csv("new_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt",sep='\t')
read_list = list(set(new_ins_reads[3]))
for ins in list(set(ins_region_reads[0])):
    sub_ins_genestat = ins_genestat[ins_genestat['newid_xl'] == ins].drop_duplicates()
    sample_gene = '|'.join(list(sub_ins_genestat['gene']))
    sample_genetype = '|'.join(list(sub_ins_genestat['genetype']))
    sample_shared = len(list(set(sub_ins_genestat['sample'])))
    if ins in list(exon_ins_region_reads[0]):
        sub_exon_ins_region_reads = exon_ins_region_reads[exon_ins_region_reads[0] == ins]
        sub_exon_ins_region_reads = sub_exon_ins_region_reads[sub_exon_ins_region_reads[3].isin(list(set(
                                    ins_region_reads[ins_region_reads[0]==ins][3])))]
        if not sub_exon_ins_region_reads.empty:
            sub_exon_ins_region_reads_gene = list(set(sub_exon_ins_region_reads[10]))
            for genename in sub_exon_ins_region_reads_gene:
                result = Counter(sub_exon_ins_region_reads[sub_exon_ins_region_reads[10] == genename][[3,10,13]].drop_duplicates()[13])
                newresult =[]
                for r in result:
                    newresult.append('exon '+str(r)+':'+str(result[r]))
                newresult = '|'.join(newresult)
                newdf = pd.DataFrame({"chr": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[0]],
                                      "start": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[1]],
                                      "end": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[2]],
                                      "INS": [ins], 'INS_type': [ins_type_all[ins_type_all[0]==ins].values.tolist()[0][-2]],
                                      'support_read_count': [len(list(set(ins_region_reads[ins_region_reads[0]==ins][3])))],
                                      'gene': [genename],
                                      'genetype': [sub_exon_ins_region_reads[sub_exon_ins_region_reads[10] == genename][11].tolist()[0]],
                                      'exon_count': [len(list(set(sub_exon_ins_region_reads[sub_exon_ins_region_reads[10] == genename][3])))],
                                      'exon_count_detail': [newresult],
                                      'sample_shared':[sample_shared],'sample_gene':[sample_gene],'sample_genetype':[sample_genetype]})
        else:
            newdf = pd.DataFrame({"chr": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[0]],
                                  "start": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[1]],
                                  "end": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[2]],
                                  "INS": [ins],
                                  'INS_type': [ins_type_all[ins_type_all[0] == ins].values.tolist()[0][-2]],
                                  'support_read_count': [
                                      len(list(set(ins_region_reads[ins_region_reads[0] == ins][3])))],
                                  'gene': ['-'],
                                  'genetype': ['-'],
                                  'exon_count': ['-'],
                                  'exon_count_detail': ['-'],
                                  'sample_shared': [sample_shared], 'sample_gene': [sample_gene],
                                  'sample_genetype': [sample_genetype]})
    else:
        newdf = pd.DataFrame({"chr": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[0]],
                              "start": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[1]],
                              "end": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[2]],
                              "INS": [ins], 'INS_type': [ins_type_all[ins_type_all[0] == ins].values.tolist()[0][-2]],
                              'support_read_count': [len(list(set(ins_region_reads[ins_region_reads[0] == ins][3])))],
                              'gene': ['-'],
                              'genetype': ['-'],
                              'exon_count': ['-'],
                              'exon_count_detail': ['-'],
                              'sample_shared':[sample_shared],'sample_gene':[sample_gene],'sample_genetype':[sample_genetype]})
    ins_ont_stat = ins_ont_stat.append(newdf)

ins_ont_stat.to_csv("targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.stat.txt",sep='\t',index = 0)


samfile = pysam.AlignmentFile("sorted_A549_ONT_rna-a549_similar_INS.bam", "rb")
outbam =  pysam.AlignmentFile("a549_ont_support_RNA_reads_stat_genename-readsName_MQ0.bam", "wb", template=samfile)
allreads=samfile.fetch()
n=0
for read in allreads:
    if read.query_name in read_list:
        outbam.write(read)
        n +=1


samfile.close()
outbam.close()
os.system(samtools + ' index a549_ont_support_RNA_reads_stat_genename-readsName_MQ0.bam')

samfile = pysam.AlignmentFile("sorted_A549_ONT_rna-a549_similar_INS.bam", "rb")
outbam =  pysam.AlignmentFile("a549_ont_support_RNA_reads_stat_genename-readsName_exon.bam", "wb", template=samfile)
allreads=samfile.fetch()
n=0
for read in allreads:
    if read.query_name in list(set(exon_ins_region_reads[3])):
        outbam.write(read)
        n +=1


samfile.close()
outbam.close()
os.system(samtools + ' index a549_ont_support_RNA_reads_stat_genename-readsName_exon.bam')

## exon length and sequence

hg38 = '/data/fs01/wangzf/nanopore/ref/mainChr_hg38_UCSC.fa'

target_INS_df = pd.read_csv("targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bed",sep='\t',header=None,
                            names=['refname','rna_satrt','rna_end','rna_read','info','strand','refname1','ins_start','ins_end'])
target_INS_exon_df = pd.read_csv("exon_targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed",sep='\t',header=None,
                                 names=['refname','rna_satrt','rna_end','rna_read','info','strand',
                                        'refname1','ins_start','ins_end','exon','gene','genetype','ENST','exon_num'])
INS_rna_df = pd.read_csv('sorted_A549_ONT_rna-a549_similar_INS.bed',sep='\t',header=None,
                         names=['refname','rna_satrt','rna_end','rna_read','info','strand'])
gene_bed = pd.read_csv('/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_noanno.gtf', sep='\t', header=None)
gene_bed['gene'] = gene_bed[8].str.split('gene_name "').str[1]
gene_bed['gene'] = gene_bed['gene'].str.split('"').str[0]
gene_bed['gene_id'] = gene_bed[8].str.split('gene_id "').str[1]
gene_bed['gene_id'] = gene_bed['gene_id'].str.split('"').str[0]
gene_bed['transcript_id'] = gene_bed[8].str.split('transcript_id "').str[1]
gene_bed['transcript_id'] = gene_bed['transcript_id'].str.split('"').str[0]
gene_bed['exon_number'] = gene_bed[8].str.split('exon_number ').str[1]
gene_bed['exon_number'] = gene_bed['exon_number'].str.split(';').str[0]
exon_bed = gene_bed[gene_bed[2]=='exon']
exon_bed1 = exon_bed.sort_values(by=['gene',3,4])
exon_bed1 = exon_bed1[[0,3,4,6,'gene','transcript_id','exon_number']]
exon_bed1.columns=['chr','start','end','strand','gene','transcript_id','exon_number']
# 需对gene_bed去重，合并完全相同的exon信息
exon_bed1['exon'] = exon_bed1['chr'] +'_' + exon_bed1['start'].map(str) + '_' + exon_bed1['end'].map(str)
exon_list = list(set(exon_bed1['exon']))
new_exon_bed1=pd.DataFrame({'chr':[],'start':[],'end':[],'strand':[],
                            'gene':[],'transcript_id':[],'exon_number':[]})
n = 0
for exon in exon_list:
    sub_exon = exon_bed1[exon_bed1['exon']==exon]
    if sub_exon.shape[0]>1:
        sub_exon1 = sub_exon.copy()
        sub_exon = sub_exon.iloc[0]
        sub_exon.loc['transcript_id'] = '|'.join(sub_exon1['transcript_id'])
        sub_exon.loc['exon_number'] = '|'.join(sub_exon1['exon_number'])
    del sub_exon['exon']
    new_exon_bed1 = new_exon_bed1.append(sub_exon)
    n +=1
    print(n)


new_exon_bed1.to_csv("/data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation_noanno.exon_remove_dup.bed",sep='\t',index=0)

out = open("A549_ONT_sniffles_smillarINS_exon_info_20210813.txt", 'w')
out1 = open("A549_ONT_sniffles_smillarINS_exon_info_exon-length_20210813.txt", 'w')
out2 = open("A549_ONT_sniffles_smillarINS_exon_info_exon-distance_20210813.txt", 'w')
out.write("#INS id\n")
out.write("#exon info: INS/hg38 chr;start of exon;end of exon;(gene;ENST-exon_number);the count of supporting reads\n")
out.write("#sequence of exon\n")
out2.write("SVid\tSV_ins\texon_type\texon_info\tdistance\n")

for refname in list(set(target_INS_df['refname'])):
    sample = refname.split("_")[0]
    hg38_loc = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
    if not os.path.exists(hg38_loc):
        hg38_loc = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
    ins_site_hg38 = open(hg38_loc).readline().strip().split()[0:3]
    out.write("\n")
    sub_target_INS_df =target_INS_df[target_INS_df['refname']==refname]
    sub_target_INS_df = sub_target_INS_df[sub_target_INS_df['info']!=0]
    sub_target_INS_df.to_csv(refname + "_all_ONT_targetINS_exon.bed", sep='\t', header=0, index=0)
    os.system(bedtools + ' merge -i {file} -d 10 -c 4,4 -o count,collapse > merged_{file}'.format(
                                        file=refname + "_all_ONT_targetINS_exon.bed"))
    merge_target_INS_exon_info = open('merged_'+refname + "_all_ONT_targetINS_exon.bed").readlines()
    if refname in list(target_INS_exon_df['refname']):
        sub_target_INS_exon_df = target_INS_exon_df[target_INS_exon_df['refname'] == refname]
        out.write('-------- '+refname + " - with adjacent gene (%s) --------\n" % (';'.join(list(set(sub_target_INS_exon_df['gene'])))))
        out1.write('-------- ' + refname + " - with adjacent gene (%s) --------\n" % (
            ';'.join(list(set(sub_target_INS_exon_df['gene'])))))
    for info in merge_target_INS_exon_info:
        info = info.strip()
        out.write(refname+"-INS\n")
        out.write(';'.join(info.split("\t"))+"\n")
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
        newinfo = info.split("\t")
        newinfo[0] = '0_segment0'
        os.system('echo "%s" > test.bed' % '\t'.join(newinfo))
        seq = os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % ins_extend).read().strip().split("\t")[-1]
        out.write(seq+"\n")
        out1.write(refname + "\tINS\t%s\n" % str(len(seq)))
    if refname in list(target_INS_exon_df['refname']):
        sub_target_INS_exon_df = target_INS_exon_df[target_INS_exon_df['refname']==refname]
        sub_target_INS_exon_df['ID'] = sub_target_INS_exon_df['ENST'] + '-' + sub_target_INS_exon_df['exon_num'].map(str)
        sub_target_INS_exon_df1 = sub_target_INS_exon_df[['exon','gene','genetype','ENST','exon_num']].drop_duplicates()
        sub_target_INS_exon_df1['ID'] = sub_target_INS_exon_df1['ENST']+'-'+sub_target_INS_exon_df1['exon_num'].map(str)
        #sub_exon_bed = exon_bed[(exon_bed['gene'].isin(list(set(sub_target_INS_exon_df1['gene'])))) & (gene_bed[2]=='exon') &
        #                        (exon_bed['transcript_id'].isin(list(set(sub_target_INS_exon_df1['ENST']))))]
        #sub_exon_bed['ID'] = sub_exon_bed['transcript_id']+'-'+sub_exon_bed['exon_number'].map(str)
        #sub_exon_bed = sub_exon_bed[sub_exon_bed['ID'].isin(sub_target_INS_exon_df1['ID'])]
        #sub_exon_bed1 = sub_exon_bed[[0, 3, 4, 'gene', 'ID']]
        sub_gene_bed = gene_bed[(gene_bed['gene'].isin(list(set(sub_target_INS_exon_df1['gene'])))) & (gene_bed[2]=='exon') &
                                (gene_bed['transcript_id'].isin(list(set(sub_target_INS_exon_df1['ENST']))))]
        sub_gene_bed['ID'] = sub_gene_bed['transcript_id']+'-'+sub_gene_bed['exon_number'].map(str)
        sub_gene_bed = sub_gene_bed[sub_gene_bed['ID'].isin(sub_target_INS_exon_df1['ID'])]
        sub_gene_bed1 = sub_gene_bed[[0, 3, 4, 'gene', 'ID']]
        exon_list = []
        for i in range(0, sub_gene_bed1.shape[0]):
            info = sub_gene_bed1.iloc[i].map(str).tolist()
            if ';'.join(info) not in exon_list:
                out.write(refname + "-hg38_gene_exon (%s)\n" % (info[-1].split('-')[-1]))
                out.write(';'.join(info) + ';' +
                          str(len(list(set(sub_target_INS_exon_df[sub_target_INS_exon_df['ID']==info[-1]]['rna_read'])))) + "\n")
                os.system('echo "%s" > test.bed' % '\t'.join(info))
                seq = \
                    os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % hg38).read().strip().split("\t")[-1]
                out.write(seq + "\n")
                for idx in merge_target_INS_exon_info:
                    ins_site_hg38[1:3] = [int(ins_site_hg38[1]), int(ins_site_hg38[2])]
                    hg38_exon_region = info[0:3]
                    hg38_exon_region[1:3] = [int(hg38_exon_region[1]), int(hg38_exon_region[2])]
                    if (hg38_exon_region[1]>ins_site_hg38[1]) and (hg38_exon_region[2]>ins_site_hg38[2]):
                        distance = hg38_exon_region[1]-ins_site_hg38[2]
                    elif (hg38_exon_region[1]<ins_site_hg38[1]) and (hg38_exon_region[2]<ins_site_hg38[2]):
                        distance=hg38_exon_region[2]-ins_site_hg38[1]
                    else:
                        distance = 0
                    outinfo = [refname, ';'.join(idx.strip().split('\t')), 'hg38_exon', ';'.join(info),str(distance)]
                    out2.write('\t'.join(outinfo)+'\n')
                exon_list.append(';'.join(info))
    sub_INS_rna_df = INS_rna_df[(INS_rna_df['refname']==refname) &
                                (INS_rna_df['rna_read'].isin(list(set(sub_target_INS_df['rna_read']))))]
    sub_INS_rna_df = sub_INS_rna_df[sub_INS_rna_df['info']>=1]
    sub_INS_rna_df.to_csv(refname+"_all_ONT_exon.bed", sep='\t',header=0,index=0)
    os.system(bedtools  + ' merge -i {file} -c 4 -o collapse > merged_{file}'.format(
                            file=refname+"_all_ONT_exon.bed"))
    merge_INS_exon_info = os.popen(bedtools +' intersect -a merged_%s -b merged_%s -v' % (refname+"_all_ONT_exon.bed",
                                    refname + "_all_ONT_targetINS_exon.bed")).read().strip().split('\n')
    for info in merge_INS_exon_info:
        info=info.strip()
        if len(list(set(info.split("\t")[-1].split(','))&set(sub_target_INS_df['rna_read']))) >=3:
            out.write(refname + "\n")
            exon_info_line = ';'.join(info.split("\t")[0:3])+';'+ \
                             str(len(list(set(info.split("\t")[-1].split(','))&set(sub_target_INS_df['rna_read']))))
            out.write(exon_info_line+ "\n")
            sample = refname.split("_")[0]
            ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
            if not os.path.exists(ins_extend):
                ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
            newinfo = info.split("\t")
            newinfo[0] = '0_segment0'
            os.system('echo "%s" > test.bed' % '\t'.join(newinfo))
            seq = \
            os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % ins_extend).read().strip().split("\t")[-1]
            out.write(seq + "\n")
            for idx in merge_target_INS_exon_info:
                ins_exon_region = idx.split('\t')[0:3]
                ins_exon_region[1:3] = [int(ins_exon_region[1]), int(ins_exon_region[2])]
                read_exon_region = info.split('\t')[0:3]
                read_exon_region[1:3] = [int(read_exon_region[1]), int(read_exon_region[2])]
                if (read_exon_region[1]>ins_exon_region[1]) and (read_exon_region[2]>ins_exon_region[2]):
                    distance = read_exon_region[1]-ins_exon_region[2]
                elif (read_exon_region[1]<ins_exon_region[1]) and (read_exon_region[2]<ins_exon_region[2]):
                    distance=read_exon_region[2]-ins_exon_region[1]
                else:
                    distance = 0
                outinfo = [refname, ';'.join(idx.strip().split('\t')), 'read_exon', exon_info_line,str(distance)]
                out2.write('\t'.join(outinfo)+'\n')
    out.write("\n")
    out2.write('\n')

out.close()
out1.close()
out2.close()
os.system("rm test.bed *_all_ONT_targetINS_exon.bed *_all_ONT_exon.bed")

#TBCE
out = open("A549_ONT_sniffles_smillarINS_exon_info_20210927-TBCE.txt", 'w')
out.write("#INS id\n")
out.write("#exon info: INS/hg38 chr;start of exon;end of exon;(gene;ENST-exon_number);the count of supporting reads\n")
out.write("#sequence of exon\n")

for refname in ['Lung94C_chr1_235409152_235410496_580']:
    sample = refname.split("_")[0]
    hg38_loc = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
    if not os.path.exists(hg38_loc):
        hg38_loc = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
    ins_site_hg38 = open(hg38_loc).readline().strip().split()[0:3]
    sub_target_INS_df =target_INS_df[target_INS_df['refname']==refname]
    sub_target_INS_df = sub_target_INS_df[sub_target_INS_df['info']!=0]
    sub_target_INS_df.to_csv(refname + "_all_ONT_targetINS_exon.bed", sep='\t', header=0, index=0)
    os.system(bedtools + ' merge -i {file} -d 10 -c 4,4 -o count,collapse > merged_{file}'.format(
                                        file=refname + "_all_ONT_targetINS_exon.bed"))
    merge_target_INS_exon_info = open('merged_'+refname + "_all_ONT_targetINS_exon.bed").readlines()
    out.write('-------- '+refname + " - with adjacent gene (TBCE) --------\n" )
    for info in merge_target_INS_exon_info:
        info = info.strip()
        out.write(refname+"-INS\n")
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
        newinfo = info.split("\t")
        newinfo[0] = '0_segment0'
        os.system('echo "%s" > test.bed' % '\t'.join(newinfo))
        seq = os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % ins_extend).read().strip().split("\t")[-1]
        out.write(seq+"\n")
        sub_target_INS_exon_df = target_INS_exon_df[target_INS_exon_df['refname']==refname]
        sub_target_INS_exon_df['ID'] = sub_target_INS_exon_df['ENST'] + '-' + sub_target_INS_exon_df['exon_num'].map(str)
        exon_list = {}
        reads = list(set(info.split('\t')[-1].split(',')))
        for read in reads:
            if read in list(sub_target_INS_exon_df['rna_read']):
                exons = sub_target_INS_exon_df[(sub_target_INS_exon_df['gene']=='TBCE') &
                                              (sub_target_INS_exon_df['rna_read']==read)]['ID'].tolist()
                for e in exons:
                    if e in list(exon_list.keys()):
                        exon_list[e]+=1
                    else:
                        exon_list[e]=1
        if exon_list:
            newline='\t'.join(info.split('\t')[0:4])+'\n'
            out.write(newline)
            newline=''
            for e in list(exon_list.keys()):
                newline += '\t'+e+':'+str(exon_list[e])
            newline +='\n'
            out.write(newline)
        out.write('\n')

out.close()


# 不涉及邻近基因外显子的新外显子信息：
for refname in ['Lung121C_chr1_29055289_29056084_159', 'Lung121C_chr1_92701576_92702371_334',
                'Lung94C_chr1_235409152_235410496_580']:
    sample = refname.split("_")[0]
    sub_target_INS_df =target_INS_df[target_INS_df['refname']==refname]
    sub_target_INS_df = sub_target_INS_df[sub_target_INS_df['info']!=0]
    sub_INS_rna_df = INS_rna_df[(INS_rna_df['refname'] == refname) &
                                (INS_rna_df['rna_read'].isin(list(set(sub_target_INS_df['rna_read']))))]
    sub_INS_rna_df = sub_INS_rna_df[sub_INS_rna_df['info'] >= 1]
    RNAMES = ' -e '.join(list(set(sub_INS_rna_df['rna_read'])))
    os.system(
        "({samtools} view -H {bam}; {samtools} view {bam} | grep -e {names}) | {samtools} view -b - {refname} > {out}".format(
            samtools=samtools, bam='sorted_A549_ONT_rna-a549_similar_INS.bam', names=RNAMES, refname=refname, out=refname+'_new_exon.bam'))
    os.system(stringtie_2_1_6 + ' -p 20 -o {refname}_new_exon.gtf {refname}_new_exon.bam'.format(refname=refname))
    os.system("cat %s_new_exon.gtf | grep exon | cut -f1,4,5,9 | cut -f1 -d ';' | \
                awk '{print $1, $2, $3, $5}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' > %s_new_exon.bed" % (refname, refname))



# 插入位点上下游最邻近的外显子
ins_bk = '/data/fs09/wangzf/nanopore/xl/A549-new/INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.500bp.bed'
hg38_exon_bed = '/data/fs08/wangzf/nanopore/xl/ref/ref/sorted_gencode.v33.annotation_exon_20200312.bed'
hg38 = '/data/fs01/wangzf/nanopore/ref/mainChr_hg38_UCSC.fa'
os.system(bedtools +' closest -a %s -b %s -D ref > closest_exon_INS_gene_INS_stat_unmapped_disc_split_split_20210705.500bp.bed' % (
            ins_bk, hg38_exon_bed))
df = pd.read_csv('closest_exon_INS_gene_INS_stat_unmapped_disc_split_split_20210705.500bp.bed',sep='\t', header=None)
df1 = df[df[13].isin(['EVI5','EPB41','TBCE'])]
df1.index = list(range(0,df1.shape[0]))
df1 = df1[[1,2,4,7,8,9,11,12,13,14,15]].drop_duplicates()
out = open('INS_closest_exon_A549_EVI5-EPB41-TBCE.txt','w')
for i in list(range(0, df1.shape[0])):
    line = '\t'.join(df1.iloc[i].map(str).tolist()[2:12])
    ense = df1.iloc[i].tolist()[-2]
    ins_site_hg38 = [int(df1.iloc[i].tolist()[0]), int(df1.iloc[i].tolist()[1])]
    info = os.popen('grep "%s\>" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf'%ense).read().strip()
    enst = info.split('transcript_id "')[1].split('"')[0]
    exon_number = info.split('exon_number ')[1].split(';')[0]
    exon1 = os.popen('grep "%s\>" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf | \
                      grep "exon_number %s" '%(enst, str(int(exon_number)-1))).read().strip()
    if exon1:
        exon_info = exon1.split('\t')
        exon_region = [int(exon_info[3]), int(exon_info[4])]
        if (exon_region[0] > ins_site_hg38[0]) and (exon_region[1] > ins_site_hg38[1]):
            distance = exon_region[0] - ins_site_hg38[1]
        elif (exon_region[1] < ins_site_hg38[0]) and (exon_region[1] < ins_site_hg38[1]):
            distance = exon_region[1] - ins_site_hg38[0]
        else:
            distance = 0
        newinfo = [df1.iloc[i].tolist()[2], exon_info[0],exon_info[3],exon_info[4],
                   exon_info[-1].split('gene_id "')[1].split('"')[0],
                   df1.iloc[i].tolist()[7],df1.iloc[i].tolist()[8],
                   exon_info[-1].split('exon_id "')[1].split('"')[0],
                   str(distance), 'former_of_ins_closest_exon']
        os.system('echo "%s" > test.bed' % ('\t'.join(newinfo[1:4])))
        fa = os.popen(bedtools +' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
        newinfo.append(fa)
        out.write('\t'.join(newinfo)+'\n')
    os.system('echo "%s" > test.bed' % ('\t'.join(df1.iloc[i].map(str).tolist()[3:6])))
    fa = os.popen(bedtools + ' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
    out.write(line+'\tins_closest_exon\t'+fa + '\n')
    exon1 = os.popen('grep "%s\>" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf | \
                      grep "exon_number %s" '%(enst, str(int(exon_number)+1))).read().strip()
    if exon1:
        exon_info = exon1.split('\t')
        exon_region = [int(exon_info[3]), int(exon_info[4])]
        if (exon_region[0] > ins_site_hg38[0]) and (exon_region[1] > ins_site_hg38[1]):
            distance = exon_region[0] - ins_site_hg38[1]
        elif (exon_region[1] < ins_site_hg38[0]) and (exon_region[1] < ins_site_hg38[1]):
            distance = exon_region[1] - ins_site_hg38[0]
        else:
            distance = 0
        newinfo = [df1.iloc[i].tolist()[2], exon_info[0],exon_info[3],exon_info[4],
                   exon_info[-1].split('gene_id "')[1].split('"')[0],
                   df1.iloc[i].tolist()[7],df1.iloc[i].tolist()[8],
                   exon_info[-1].split('exon_id "')[1].split('"')[0],
                   str(distance), 'latter_of_ins_closest_exon']
        os.system('echo "%s" > test.bed' % ('\t'.join(newinfo[1:4])))
        fa = os.popen(bedtools +' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
        newinfo.append(fa)
        out.write('\t'.join(newinfo)+'\n')


out.close()

ins_exon= open('A549_ONT_sniffles_smillarINS_exon_info_20210813.txt')
out = open('A549_ONT_sniffles_smillarINS_exon_info_20210813-nearExon.txt','w')
os.mkdir('./EPB41_EVI5_TBCE_gtf')
line = ins_exon.readline()
while line:
    if 'hg38_gene_exon' in line:
        sv = line.split('-')[0]
        line1 =ins_exon.readline()
        info = line1.split(';')
        gene = info[3]
        enst = info[4].split('-')[0]
        if gene in ['EPB41','EVI5','TBCE']:
            os.system('grep "%s" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf | \
                       grep "%s" > ./EPB41_EVI5_TBCE_gtf/%s-%s.gtf' %
                      (gene, enst, gene, enst))
        exon_number = info[4].split('-')[1]
        exon1 = os.popen('grep "%s\>" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf | \
                          grep "exon_number %s" ' % (enst, str(int(exon_number) - 1))).read().strip()
        if exon1:
            exon_info = exon1.split('\t')
            exon_region = [int(exon_info[3]), int(exon_info[4])]
            newinfo = [sv, exon_info[0], exon_info[3], exon_info[4],
                       gene, enst,'exon number ' +str(int(exon_number) - 1),
                       'former_of_ins_link_exon']
            os.system('echo "%s" > test.bed' % ('\t'.join(newinfo[1:4])))
            fa = os.popen(bedtools + ' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
            newinfo.append(fa)
            out.write('\t'.join(newinfo) + '\n')
        os.system('echo "%s" > test.bed' % ('\t'.join(info[0:3])))
        fa = os.popen(bedtools + ' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
        out.write(line + '\tins_closest_exon\t' + fa + '\n')
        exon1 = os.popen('grep "%s\>" /data/fs08/wangzf/nanopore/xl/ref/ref/gencode.v33.annotation.gtf | \
                          grep "exon_number %s" ' % (enst, str(int(exon_number) + 1))).read().strip()
        if exon1:
            exon_info = exon1.split('\t')
            exon_region = [int(exon_info[3]), int(exon_info[4])]
            newinfo = [sv, exon_info[0], exon_info[3], exon_info[4],
                       gene, enst,'exon number ' +str(int(exon_number) + 1),
                       'latter_of_ins_link_exon']
            os.system('echo "%s" > test.bed' % ('\t'.join(newinfo[1:4])))
            fa = os.popen(bedtools + ' getfasta -fi %s -bed test.bed' % hg38).read().strip().split('\n')[-1]
            newinfo.append(fa)
            out.write('\t'.join(newinfo) + '\n')
        out.write('\n')
    line = ins_exon.readline()




ins_exon.close()
out.close()

os.system('tar -zcvf ./EPB41_EVI5_TBCE_gtf.tar.gz EPB41_EVI5_TBCE_gtf/')
os.system('aws s3 cp EPB41_EVI5_TBCE_gtf.tar.gz s3://for_igv --endpoint-url http://10.100.116.91:9020')

# TBCE, EPB41, EVI5
refnames = {'EPB41':['Lung121C_chr1_29055289_29056084_159'],
            #'EVI5':['Lung121C_chr1_92701576_92702371_334'],
            'TBCE':['Lung94C_chr1_235409152_235410496_580']}
target_INS_df = pd.read_csv("targetRegion_sorted_A549_ONT_rna-a549_similar_INS.bed",sep='\t',header=None,
                            names=['refname','rna_satrt','rna_end','rna_read','info','strand','refname1','ins_start','ins_end'])
target_INS_exon_df = pd.read_csv("exon_targetRegion_sorted_A549_ONT_rna-a549_similar_INS_reads.bed",sep='\t',header=None,
                                 names=['refname','rna_satrt','rna_end','rna_read','info','strand',
                                        'refname1','ins_start','ins_end','exon','gene','genetype','ENST','exon_num'])
os.system(samtools +' index -@ 20 sorted_A549_ONT_rna-a549_similar_INS.bam')
for gene in list(refnames.keys()):
    for refname in refnames[gene]:
        sample = refname.split("_")[0]
        hg38_loc = os.path.join(ins_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
        if not os.path.exists(hg38_loc):
            hg38_loc = os.path.join(dup_locdir, sample, refname, refname + '_SRM_INS_hg38_coordinate.bed')
        ins_site_hg38 = open(hg38_loc).readline().strip().split()[0:3]
        sub_target_INS_df =target_INS_df[target_INS_df['refname']==refname]
        sub_target_INS_df = sub_target_INS_df[sub_target_INS_df['info']!=0]
        sub_target_INS_exon_df = target_INS_exon_df[target_INS_exon_df['refname']==refname]
        reads = list(set(sub_target_INS_exon_df[sub_target_INS_exon_df['rna_read'].isin(sub_target_INS_df['rna_read'])]['rna_read']))
        new_sub_target_INS_exon_df = sub_target_INS_exon_df[sub_target_INS_exon_df['rna_read'].isin(reads)]
        enst_list = list(set(new_sub_target_INS_exon_df['ENST']))
        c_enst_list = Counter(list(new_sub_target_INS_exon_df['ENST']))
        # reads 提取bam
        for enst in enst_list:
            if c_enst_list[enst] >= 3:
                newreads = list(set(new_sub_target_INS_exon_df[new_sub_target_INS_exon_df['ENST']==enst]['rna_read']))
                ReadName = ' -e '.join(newreads)
                os.system(
                    "({samtools} view -H {bam}; {samtools} view {bam} {refname} | grep -e {names}) | {samtools} view -bS - > {out}".format(
                        samtools=samtools, bam='sorted_A549_ONT_rna-a549_similar_INS.bam', names=ReadName, refname=refname,
                        out= './target_A549_new/' + '_'.join([gene, refname, enst]) + '_new_exon_reads.bam'))
                os.system(samtools + ' index -@ 10 ./target_A549_new/' + '_'.join([gene, refname, enst]) + '_new_exon_reads.bam')
                os.system('cp %s* ./target_A549_new/new_exon_3gene/' % ('./target_A549_new/' + '_'.join([gene, refname,enst]) + '_new_exon_reads.bam'))
        # 覆盖度最高的INS区域与联动的exon
        # INS区域
        sub_target_INS_df1 = sub_target_INS_df[sub_target_INS_df['rna_read'].isin(reads)]
        subdir = gene+'_'+refname
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        file_list = []
        for read in list(set(sub_target_INS_df1['rna_read'])):
            subdf = sub_target_INS_df1[sub_target_INS_df1['rna_read']==read]
            file = '-'.join([gene,refname,read])+'.bed'
            subdf.to_csv(os.path.join('.',subdir,file),sep='\t',header=None,index=None)
            file_list.append(os.path.join('.',subdir,file))
        os.system(bedtools + " multiinter -i %s | awk -F '\t' '{if ($4>=2) print $0}' - | %s merge -i - -d 10 | \
                  awk -F '\t' '{if ($3-$2>30) print $1\"\t\"$2\"\t\"$3\"\t\"$3-$2}' - > %s" %
                  (' '.join(file_list), bedtools, '_'.join([gene,refname,'merged_cov_region_INSexon.bed'])))
        # exon区域，覆盖度大于2
        sub_target_INS_exon_df1 = sub_target_INS_exon_df[sub_target_INS_exon_df['rna_read'].isin(reads)]
        file_list = []
        for read in list(set(sub_target_INS_exon_df1['rna_read'])):
            subdf = sub_target_INS_exon_df1[sub_target_INS_exon_df1['rna_read'] == read]
            file = '-'.join([gene, refname, read]) + '_exon.bed'
            subdf.to_csv(os.path.join('.', subdir, file), sep='\t', header=None, index=None)
            file_list.append(os.path.join('.', subdir, file))
        os.system(bedtools + " multiinter -i %s | awk -F '\t' '{if ($4>=2) print $0}' - | %s merge -i - -d 10 | \
                          awk -F '\t' '{if ($3-$2>30) print $1\"\t\"$2\"\t\"$3\"\t\"$3-$2}' - > %s" %
                  (' '.join(file_list), bedtools, '_'.join([gene, refname, 'merged_cov_region_exon.bed'])))
        os.system('rm -r '+ subdir)
        # 提取序列
        outfile = os.path.join('target_A549_new/new_exon_3gene', '_'.join([gene,refname,'_seqInfo_20211012.txt']))
        os.system('grep "%s" a549_similar_somaticINS_ins_v2n.bed | awk -F \'\t\' \'{print $0"\tINS region"}\' > %s' % (refname, outfile))
        out = open(outfile,'a')
        # INS序列
        ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
        if not os.path.exists(ins_extend):
            ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
        os.system('grep "%s" a549_similar_somaticINS_ins_v2n.bed |awk -F \'\t\' \'{print "0_segment0\t"$2"\t"$3}\' - > test.bed' % refname)
        seq = os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % ins_extend).read().strip().split("\t")[-1]
        out.write('\t'.join([gene, refname, 'INS.sequence'])+'\n')
        out.write(seq+'\n')
        out.close()
        os.system('cat %s >> %s' % ('_'.join([gene, refname, 'merged_cov_region_INSexon.bed']), outfile))
        # exon序列
        out = open(outfile, 'a')
        out.write('\n')
        shared_exon_seq = os.popen("awk -F '\t' '{print \"0_segment0\t\"$2\"\t\"$3}' %s |" % ('_'.join([gene, refname, 'merged_cov_region_exon.bed'])) +
                                   bedtools + ' getfasta -fi %s -bed - -tab' % (ins_extend)).read().strip().split('\n')
        for seq in shared_exon_seq:
            out.write(refname + '\t' + seq.split(':')[1].split('\t')[0] + ' ---- exon \n')
            out.write(seq.split(':')[1].split('\t')[1]+'\n')
        out.close()
        os.system('rm test.bed')



## NGS
a549_NGSrna_fq1_1 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate01_1.fastq.gz'
a549_NGSrna_fq1_2 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate01_2.fastq.gz'
a549_NGSrna_fq2_1 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate02_1.fastq.gz'
a549_NGSrna_fq2_2 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate02_2.fastq.gz'
a549_NGSrna_fq3_1 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate03_1.fastq.gz'
a549_NGSrna_fq3_2 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_replicate03_2.fastq.gz'


a549_NGSrna_fq1 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_allreplicate_1.fastq.gz'
a549_NGSrna_fq2 = '/data/fs01/wangzf/nanopore/cell_line/RNA_ascp/A549_allreplicate_2.fastq.gz'
os.system('zcat %s %s %s | gzip - > %s' % (a549_NGSrna_fq1_1, a549_NGSrna_fq2_1, a549_NGSrna_fq3_1, a549_NGSrna_fq1))
os.system('zcat %s %s %s | gzip - > %s' % (a549_NGSrna_fq1_2, a549_NGSrna_fq2_2, a549_NGSrna_fq3_2, a549_NGSrna_fq2))
# hg38-NGS
os.system(hisat2_build + " %s %s -p 10" % (hg38, workdir))
sam = 'A549_NGS_rna.sam'
os.system(hisat2 + ' -p 10 -x {sample_workdir} -1 {fq1} -2 {fq2} -S {sam}'.format(
    sample_workdir=workdir, fq1=a549_NGSrna_fq1, fq2=a549_NGSrna_fq2, sam=sam))
os.system(samtools + " view -bS %s -@ 10 > %s.bam" % (sam, sam.split(".sam")[0]))
os.system(samtools + " sort -@ 10 %s.bam > sorted_%s.bam" % (sam.split(".sam")[0], sam.split(".sam")[0]))
# os.system('rm *sam')
os.system(stringtie + " sorted_%s.bam -m 100 -c 1 -o %s.gtf -p 10" % (sam.split(".sam")[0], sam.split(".sam")[0]))
os.system(gffcompare + ' -R -r %s -o A549_stringtieAssembled_comp_gencodev33 %s.gtf' % (gencode, sam.split('.sam')[0]))
tmap = pd.read_csv("A549_stringtieAssembled_comp_gencodev33.A549_NGS_rna.gtf.tmap", sep='\t')
new_tmap = tmap[tmap['ref_gene_id'].isin(candidate_gene_ID)]
new_tmap.to_csv("filter1_A549_NGS_gencodev33_compare.tmap", sep='\t', index=0)
# 过滤，只保留exon>1并且长度>200bp的transcripts
os.system("awk '($6>1 && $10>=200){print$0}' filter1_A549_NGS_gencodev33_compare.tmap \
            > filter2_by_exon_length_A549_NGS.tmap")
os.system("awk '{print $5}' filter2_by_exon_length_A549_NGS.tmap > filter2_by_exon_length_A549_NGS_transcript_ID")
os.system("grep -Ff filter2_by_exon_length_A549_NGS_transcript_ID -w A549_ONT_rna.gtf > \
            filter2_by_exon_length_A549_NGS_transcript.gtf")
