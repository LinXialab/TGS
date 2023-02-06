
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

samtools = '/biosoft/samtools-1.9/samtools'
hisat2 = '/wangzf/software/hisat2-2.1.0/hisat2'
hisat2_build = '/wangzf/software/hisat2-2.1.0/hisat2-build'
stringtie = '/wangzf/software/stringtie-2.1.4/stringtie'
stringtie_2_1_6 = '/tools/stringtie-2.1.6.Linux_x86_64/stringtie'
seqtk = '/tools/seqtk/seqtk'
ngmlr = "/software/ngmlr-0.2.7/ngmlr"
hg38 = '/nanopore/ref/mainChr_hg38_UCSC.fa'
ins_locdir = '/INS_RepeatMasker_f50bp/Assembly_result'
dup_locdir = '/Assembly_DUP/DUP_2286/Assembly_result'
gencode = '/ref/gencode.v33.annotation.gtf'

workdir = '/PRJDB9872'
# NANOPORE RNA: DRX218781: MinION sequencing of SAMD00226192(Aberrant transcript isoforms detected by full-length transcriptome sequencing as transcripts of potential neoantigens in non-small cell lung cancer)
# H1299 & PC-9 RNA # 10.100.2.3
minimap2 = '/software/minimap2/minimap2'
gffcompare = '/software/gffcompare-0.12.2.Linux_x86_64/gffcompare'

H1299_ONTrna_fq = '/PRJDB9872/DRR228519_1.fastq.gz'
PC_9_ONTrna_fq = '/PRJDB9872/DRR228530_1.fastq.gz'
cell_fq = {'H1299':H1299_ONTrna_fq, 'PC-9':PC_9_ONTrna_fq}

if not os.path.exists(os.path.join(workdir, 'ONT_NGS_RNA')):
    os.mkdir(os.path.join(workdir, 'ONT_NGS_RNA'))


os.chdir(os.path.join(workdir, 'ONT_NGS_RNA'))

# all INS fasta
ins = pd.read_csv('/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all_samples_5709_newINS.bed',
                  sep='\t',header=None)
new_fasta = open('/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all_samples_5709_newINS.fasta','w')
for refname in list(set(ins[7])):
    sample = refname.split('_')[0]
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
                'awk -F \'\t\' \'{print "%s\t"$2"\t"$3}\' %s >> all_somaticINS_v2n.bed' % (refname, locfile))
            n = 1
            for line in open(ins_extend).readlines():
                if ">" in line:
                    if n == 1:
                        newline = '>' + refname + '\n'
                        n += 1
                    else:
                        newline = '>' + refname + '-' + str(n) + '\n'
                        n += 1
                    new_fasta.write(newline)
                else:
                    new_fasta.write(line)


new_fasta.close()
os.system('sort -k1,1 -k2,2n all_somaticINS_v2n.bed > sorted_all_somaticINS_v2n.bed')
ins_fasta = '/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all_samples_5709_newINS.fasta'
## ONT INS
for cell in list(cell_fq.keys()):
    ONTrna_fq = cell_fq[cell]
    os.system(minimap2 + ' -ax splice -uf -k14 {} {} >\
                            {}_ONT_rna_INS.sam'.format(ins_fasta, ONTrna_fq, cell))
    os.system(samtools + " view -bS {}_ONT_rna_INS.sam -@ 10 > {}_ONT_rna_INS.bam".format(cell, cell))
    os.system(samtools + " sort {}_ONT_rna_INS.bam > sorted_{}_ONT_rna_INS.bam".format(cell, cell))
    os.system(samtools + ' index sorted_{}_ONT_rna_INS.bam'.format(cell))
    os.system(bedtools + ' bamtobed -split -i sorted_{}_ONT_rna_INS.bam | sort -k1,1 -k2,2n > \
                            sorted_{}_ONT_rna_INS.bed'.format(cell, cell))
    # if ONT reads mapped on INS seq overlapped INS region
    os.system(bedtools + ' intersect -a sorted_{}_ONT_rna_INS.bed \
                            -b sorted_all_somaticINS_v2n.bed -wa -wb > \
                            targetRegion_sorted_{}_ONT_rna_INS.bed'.format(cell, cell)) # read bed overlapped with INS regions
    target_INS = pd.read_csv("targetRegion_sorted_{}_ONT_rna_INS.bed".format(cell),sep='\t',header=None,
                             usecols = [0,1,2,3],names=['refname','start','end','read'])
    os.system('/biosoft/sambamba-0.7.0/sambamba slice -L sorted_all_somaticINS_v2n.bed \
                sorted_{}_ONT_rna_INS.bam > targetRegion_sorted_{}_ONT_rna_INS.bam'.format(cell, cell))
    reads = list(set(target_INS['read']))
    samfile = pysam.AlignmentFile("sorted_{}_ONT_rna_INS.bam".format(cell), "rb")
    outbam =  pysam.AlignmentFile(cell + "_ont_support_RNA_reads_stat_genename-readsName.bam", "wb", template=samfile)
    allreads=samfile.fetch()
    n=0
    for read in allreads:
        if read.query_name in reads:
            outbam.write(read)
            n +=1
    samfile.close()
    outbam.close()
    os.system(samtools + ' index %s_ont_support_RNA_reads_stat_genename-readsName.bam' % cell)
    os.system(stringtie + " sorted_{}_ONT_rna_INS.bam -m 100 -c 1 -o sorted_{}_ONT_rna_INS.gtf -p 10".format(cell, cell))
    os.system(gffcompare + ' -R -r {} -o  sorted_{}_INS_comp_gencodev33 sorted_{}_ONT_rna_INS.gtf'.format(gencode, cell, cell))


# if ONT reads mapped on INS region mapped on known exon
newgtf = '/INS_DUP_RNA_20210630/INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.gtf'
gtf_df = pd.read_csv(newgtf, sep='\t',header=None)
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
exon_gtf.to_csv("INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.exon.bed",sep='\t',index=0,header=0)

for cell in list(cell_fq.keys()):
    ins_reads = pd.read_csv("sorted_{}_ONT_rna_INS.bed".format(cell),sep='\t',header=None)
    ins_region_reads = pd.read_csv("targetRegion_sorted_{}_ONT_rna_INS.bed".format(cell),sep='\t',header=None)
    new_ins_reads = ins_reads[ins_reads[3].isin(list(set(ins_region_reads[3])))]
    new_ins_reads.to_csv('targetRegion_sorted_{}_ONT_rna_INS_reads.bed'.format(cell),sep='\t',header=0,index=0)
    os.system(bedtools +' intersect -a targetRegion_sorted_{}_ONT_rna_INS_reads.bed \
                            -b INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210715_gencodeV33.exon.bed -wa -wb >\
                             exon_targetRegion_ssorted_{}_ONT_rna_INS_reads.bed'.format(cell, cell))
    exon_ins_region_reads = pd.read_csv("exon_targetRegion_ssorted_{}_ONT_rna_INS_reads.bed".format(cell),sep='\t',header=None)
    ins_ont_stat = pd.DataFrame({"chr":[],"start":[],"end":[],"INS":[],'INS_type':[],
                                 'support_read_count':[],'gene':[],'genetype':[],
                                 'exon_count':[],'exon_count_detail':[],
                                 'sample_shared':[],'sample_gene':[],'sample_genetype':[]})
    ins_bed = pd.read_csv('/new_somatic_SV_TGS_bed_20210524/somatic_sv_new_bed/all.INS.DUP.anno.region.clean.20210524.tsv',
                          sep='\t',header=0)
    oldworkdir = '/new_somatic_SV_TGS_bed_20210524/INS_DUP_RNA_20210624'
    ins_type_all = pd.read_csv(oldworkdir + "/INS_repeat_20210624.bed", header=None, sep='\t')
    ins_genestat = pd.read_csv("/A549/ONT_NGS_RNA/new_INS_allSamples_repeat_gene_INS_stat_unmapped_disc_split_20210705.txt",sep='\t')
    read_list = list(set(new_ins_reads[3]))
    for ins in list(set(ins_region_reads[0])):
        sub_ins_genestat = ins_genestat[ins_genestat['newid_xl'] == ins].drop_duplicates()
        if not sub_ins_genestat.empty:
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
                        if ins in list(ins_type_all[0]):
                            insType = ins_type_all[ins_type_all[0] == ins].values.tolist()[0][-2]
                        else:
                            insType = '-'
                        newdf = pd.DataFrame({"chr": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[0]],
                                              "start": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[1]],
                                              "end": [ins_bed[ins_bed['newid_xl']==ins].values.tolist()[0][-2].split("_")[2]],
                                              "INS": [ins], 'INS_type': [insType],
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
                if ins in list(ins_type_all[0]):
                    insType = ins_type_all[ins_type_all[0] == ins].values.tolist()[0][-2]
                else:
                    insType = '-'
                newdf = pd.DataFrame({"chr": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[0]],
                                      "start": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[1]],
                                      "end": [ins_bed[ins_bed['newid_xl'] == ins].values.tolist()[0][-2].split("_")[2]],
                                      "INS": [ins], 'INS_type': [insType],
                                      'support_read_count': [len(list(set(ins_region_reads[ins_region_reads[0] == ins][3])))],
                                      'gene': ['-'],
                                      'genetype': ['-'],
                                      'exon_count': ['-'],
                                      'exon_count_detail': ['-'],
                                      'sample_shared':[sample_shared],'sample_gene':[sample_gene],'sample_genetype':[sample_genetype]})
            ins_ont_stat = ins_ont_stat.append(newdf)
    ins_ont_stat.to_csv("targetRegion_sorted_{}_INS_reads.stat.txt".format(cell),sep='\t',index = 0)
    samfile = pysam.AlignmentFile("sorted_{}_ONT_rna_INS.bam".format(cell), "rb")
    outbam =  pysam.AlignmentFile("{}_ont_support_RNA_reads_stat_genename-readsName_MQ0.bam".format(cell), "wb", template=samfile)
    allreads=samfile.fetch()
    n=0
    for read in allreads:
        if read.query_name in read_list:
            outbam.write(read)
            n +=1
    samfile.close()
    outbam.close()
    os.system(samtools + ' index {}_ont_support_RNA_reads_stat_genename-readsName_MQ0.bam'.format(cell))
    samfile = pysam.AlignmentFile("sorted_{}_ONT_rna_INS.bam".format(cell), "rb")
    outbam =  pysam.AlignmentFile("{}_ont_support_RNA_reads_stat_genename-readsName_exon.bam".format(cell), "wb", template=samfile)
    allreads=samfile.fetch()
    n=0
    for read in allreads:
        if read.query_name in list(set(exon_ins_region_reads[3])):
            outbam.write(read)
            n +=1
    samfile.close()
    outbam.close()
    os.system(samtools + ' index {}_ont_support_RNA_reads_stat_genename-readsName_exon.bam'.format(cell))

## exon length and sequence

hg38 = '/ref/mainChr_hg38_UCSC.fa'
gene_bed = pd.read_csv('/ref/gencode.v33.annotation_noanno.gtf', sep='\t', header=None)
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


for cell in list(cell_fq.keys()):
    target_INS_df = pd.read_csv("targetRegion_sorted_{}_ONT_rna_INS.bed".format(cell),sep='\t',header=None,
                                names=['refname','rna_satrt','rna_end','rna_read','info','strand','refname1','ins_start','ins_end'])
    target_INS_exon_df = pd.read_csv("exon_targetRegion_ssorted_{}_ONT_rna_INS_reads.bed".format(cell),sep='\t',header=None,
                                     names=['refname','rna_satrt','rna_end','rna_read','info','strand',
                                            'refname1','ins_start','ins_end','exon','gene','genetype','ENST','exon_num'])
    INS_rna_df = pd.read_csv("sorted_{}_ONT_rna_INS.bed".format(cell),sep='\t',header=None,
                             names=['refname','rna_satrt','rna_end','rna_read','info','strand'])
    out = open(cell+"_ONT_all_INS_exon_info_20211028.txt", 'w')
    out1 = open(cell+"_ONT_all_INS_exon_info_20211028_exon-length.txt", 'w')
    out2 = open(cell+"_ONT_all_INS_exon_info_20211028_exon-distance.txt", 'w')
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
        sub_target_INS_df =target_INS_df[target_INS_df['refname']==refname]
        sub_target_INS_df = sub_target_INS_df[sub_target_INS_df['info']!=0]
        if not sub_target_INS_df.empty:
            out.write("\n")
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


os.chdir('/PRJDB9872/ONT_NGS_RNA/h1299_pc-9')
os.system('cp ../H1299*txt ./')
os.system('cp ../PC-9*txt ./')
os.system("grep 'gene' H1299_ONT_all_INS_exon_info_20211028_exon-length.txt |\
            awk -F ' ' '{print $2\"\t\"$7}' | awk -F '(' '{print $1\"\t\"$2}' | awk -F ')' '{print $1}'| \
            sort -u > H1299_ONT_all_INS_exon_info_20211028_adjacentGene.txt")
os.system("grep 'gene' PC-9_ONT_all_INS_exon_info_20211028_exon-length.txt |\
            awk -F ' ' '{print $2\"\t\"$7}' | awk -F '(' '{print $1\"\t\"$2}' | \
            awk -F ')' '{print $1}'| sort -u > PC-9_ONT_all_INS_exon_info_20211028_adjacentGene.txt")

for cell in list(cell_fq.keys()):
    df = pd.read_csv(cell+'_ONT_all_INS_exon_info_20211028_adjacentGene.txt', sep='\t',header=None)
    df = df[[0,2]]
    exon_counts = open(cell+"_ONT_all_INS_exon_info_20211028.txt")
    out = open(cell+'_ONT_all_INS_exon_info_20211028_adjacentGene_filtered_supportReads.txt', 'w')
    out.write('ins\tgene\tchr\tstart_exon\tend_exon\tENST-exon_number\tsupporting_reads\n')
    line = exon_counts.readline()
    while line:
        if (';' in line) and ('ENST' in line):
            info = line.strip().split(';')
            if info[3] in list(df[2]):
                if int(info[-1])>=5:
                    newinfo = info.copy()
                    del newinfo[3]
                    newline = '\t'.join(df[df[2]==info[3]].values.tolist()[0]) + '\t' + '\t'.join(info) +'\n'
                    out.write(newline)
        line = exon_counts.readline()
    exon_counts.close()
    out.close()



gene_list = {'H1299':['MCL1','SHC4','ACTR3','TRIM4','CPSF2','HCP5','LARP6','LONP1','SFR1','SNHG20',
                      'TRA2A','TUSC3'],
             'PC-9':[#'MCL1','SFR1','MIF'
                     'ANXA2','CLNS1A','SDHAP1','SDHAP3']}

def cell_gene_bam(cell,gene):
        info_list = os.popen('grep "{}" {}_ONT_all_INS_exon_info_20211028_exon-distance.txt'.format(gene, cell)).read().strip().split('\n')
        for info in info_list:
            if info:
                refname = info.split('\t')[0]
                enst_exon = info.split('\t')[-2].split(';')[-1]
                reads = info.split('\t')[1].split(';')[-1].split(',')
                RNAMES = ' -e '.join(list(set(reads)))
                os.system(
                    "({samtools} view -H {bam}; {samtools} view {bam} {refname}| grep -e {names}) | \
                    {samtools} view -@ 10 -b - > {out}".format(
                        samtools=samtools, bam='../sorted_%s_ONT_rna_INS.bam' % cell, names=RNAMES, refname=refname,
                        out='_'.join([refname, gene, cell,enst_exon, 'support_INS_exon.bam'])))
                os.system(samtools + ' index %s' % ('_'.join([refname, gene, cell, enst_exon,'support_INS_exon.bam'])))
            else:
                print(gene+ ' ' + cell)


pools = Pool(3)
for cell in list(gene_list.keys()):
    df = pd.read_csv(cell + '_ONT_all_INS_exon_info_20211028_adjacentGene.txt', sep='\t', header=None)
    df = df[[0, 2]]
    for gene in gene_list[cell]:
        pools.apply_async(cell_gene_bam,args=(cell, gene,))


pools.close()
pools.join()
del pools


# (SFR1)
os.chdir('/PRJDB9872/ONT_NGS_RNA')
read_bed = 'sorted_PC-9_ONT_rna_INS.bed'
out = open("PC-9_ONT_sniffles_smillarINS_exon_info_20211106-SFR1.txt", 'w')
refname='Lung95C_chr10_104126288_104126379_321'
out.write(refname+' - SFR1\n')
data=open("PC-9_ONT_all_INS_exon_info_20211028.txt")
line = data.readline()
while line:
    if 'SFR1' in line:
        ins_list = {}
        ins_list[data.readline()]=[data.readline(),data.readline()]
        ins_list[data.readline()] = [data.readline(), data.readline()]
        break
    line = data.readline()


data.close()
cell='PC-9'
target_INS_exon_df = pd.read_csv("exon_targetRegion_ssorted_{}_ONT_rna_INS_reads.bed".format(cell),sep='\t',header=None,
                                     names=['refname','rna_satrt','rna_end','rna_read','info','strand',
                                            'refname1','ins_start','ins_end','exon','gene','genetype','ENST','exon_num'])
for ins in list(ins_list.keys()):
    insID = ins_list[ins][0].split('-INS')[0]
    enst_list = []
    reads = ins_list[ins][1].strip().split(';')[-1].split(',')
    for read in reads:
        subdf = target_INS_exon_df[(target_INS_exon_df['refname'] == 'Lung95C_chr10_104126288_104126379_321') & (
                                    target_INS_exon_df['rna_read'] == read)][['ENST','exon_num']]
        subdf['enstid'] = subdf['ENST'] +'-'+subdf['exon_num'].map(str)
        enst_list.extend(list(set(subdf['enstid'])))
    enstlist = list(set(enst_list))
    for line in ins_list[ins]:
        out.write(line)
    out.write('\t'.join(enstlist)+'\n')


out.close()


# h1299
os.chdir('/PRJDB9872/ONT_NGS_RNA')
read_bed = 'sorted_H1299_ONT_rna_INS.bed'
out = open("H1299_ONT_sniffles_smillarINS_exon_info_2021110.txt", 'w')
data=open("H1299_ONT_all_INS_exon_info_20211028.txt")
line = data.readline()
ins_list = {}
gene_ins = {'Lung21C_chr1_150575096_150576097_402;7455;8106;2':'MCL1',
            'Lung21C_chr1_150575096_150576097_402;9002;10520;155':'MCL1',
            'Lung21C_chr2_113948950_113948981_286;15109;15419;2':'ACTR3',
            'Lung90C_chr14_92152593_92153558_170;21293;22008;10':'CPSF2',
            'Lung100C_chr19_5707949_5708960_67;16075;16691;2':'LONP1',
            'Lung95C_chr10_104126288_104126379_321;12090;12951;36':'SFR1'
            }
while line:
    for gene in list(gene_ins.keys()):
        if gene in line:
            ins = gene.split(';')[0]+'-INS'
            ins_list[gene]=[ins,line]
    line = data.readline()


data.close()
cell='H1299'
target_INS_exon_df = pd.read_csv("exon_targetRegion_ssorted_{}_ONT_rna_INS_reads.bed".format(cell),sep='\t',header=None,
                                     names=['refname','rna_satrt','rna_end','rna_read','info','strand',
                                            'refname1','ins_start','ins_end','exon','gene','genetype','ENST','exon_num'])
for ins in list(ins_list.keys()):
    insID = ins_list[ins][0].split('-INS')[0]
    out.write(gene_ins[ins] +'\n')
    enst_list = []
    reads = ins_list[ins][1].strip().split(';')[-1].split(',')
    for read in reads:
        subdf = target_INS_exon_df[(target_INS_exon_df['refname'] ==insID) & (
                                    target_INS_exon_df['rna_read'] == read)][['ENST','exon_num']]
        if not subdf.empty:
            subdf['enstid'] = subdf['ENST'] +'-'+subdf['exon_num'].map(str)
            enst_list.extend(list(set(subdf['enstid'])))
    enstlist = list(set(enst_list))
    out.write(ins_list[ins][0]+'\n')
    out.write('\t'.join(enstlist)+'\n')
    out.write('\n')


out.close()

os.chdir('/xl/PRJDB9872/new_INS_seq_wh_20211116')
for filename in os.listdir('./'):
    if ('seq' not in filename) and ('.bed' not in filename):
        file = open(filename)
        out = open('seq_'+filename,'w')
        for line in file.readlines():
            if line != '\n':
                out.write(line)
                info = line.strip().split(';')
                sample = info[1].split('_')[0]
                refname = info[1]
                ins_extend = os.path.join(ins_locdir, sample, refname, refname + '_SRM.fasta')
                if not os.path.exists(ins_extend):
                    ins_extend = os.path.join(dup_locdir, sample, refname, refname + '_SRM.fasta')
                newinfo = info.copy()
                newinfo[1] = '0_segment0'
                os.system('echo "%s" > test.bed' % '\t'.join(newinfo[1:4]))
                seq = \
                    os.popen(bedtools + ' getfasta -fi %s -bed test.bed -tab' % ins_extend).read().strip().split("\t")[-1]
                out.write(seq + "\n")
        out.close()





















