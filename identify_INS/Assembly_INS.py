# __author__ = 'ztf-20200918'
# !/usr/bin/python
# -*- coding:utf-8 -*-

import os
import pysam
from multiprocessing import Pool
import time

samtools = "/app/samtools-1.9/samtools"
bedtools = "/tools/bedtools2/bin/bedtools"
hisat2 = "/tools/hisat2-2.1.0"
bwa = "/Tools/bin/bwa"
shasta = '/Tools/shasta-Linux-0.5.1'
minimap2 = "/Mapping/minimap2-master/minimap2"
mp_helen = "/Assembly/MarginPolish_Helen_docker.sh"
seqtk = "/tools/seqtk/seqtk"
racon = "/tools/racon/build/bin/racon"
pilon = "/tools/pilon-1.23.jar"

hg38 = '/human_hg38/GRCh38/hg38_UCSC/mainChr_hg38_UCSC.fa'
ont_bam_dir = "/TGS/Nanopore/mapping-ngmlr"
atac_bam_dir = "/ATAC-seq/bam"
rna_bam_dir = "/TGS/RNA-seq/pipeline"
ngs_bam_dir = "/WGS/mapping"
ngs_fq_dir = "/DataStore/tgs/lung_cancer/WGS/cleandata"
atac_fq_dir = "/DataStore/tgs/lung_cancer/ATAC-seq/cleandata"
rna_fq_dir = "/DataStore/tgs/lung_cancer/RNA-seq/cleandata"
ins_region_dir = "/INS/somatic_INS_select"
ins_sniffles_dir = "/INS/breakpoint20200723"

########################################################################################
# samples
########################################################################################
sample_list = ['Lung01', 'Lung02', 'Lung05', 'Lung06', 'Lung08']

work_dir = "/Assembly"
ngs_pre_dir = os.path.join(work_dir, "NGS")
asmb_dir = os.path.join(work_dir, "INS")


def after_align(sam):
    prefix = sam.split('.sam')[0]
    bam = "%s.bam" % prefix
    bam_sort = "%s.sorted.bam" % prefix
    os.system("{samtools} view -F 4 -q 30 -S -b {sam} -o {bam}".format(samtools=samtools, sam=sam, bam=bam))
    os.system("{samtools} sort {bam} -o {bams}".format(samtools=samtools, bam=bam, bams=bam_sort))
    os.system('{samtools} index {bams}'.format(bams=bam_sort, samtools=samtools))
    os.system("rm %s" % bam)
    os.system("rm %s" % sam)


########################################################################################
# extract variant ont reads from bam
########################################################################################
def variant_ont_qnames(svid_x, sv_file_x, vr_txt_x, ont_bam_x, region_x, ont_bam_region_x, vr_bam_x, vr_fq_x):
    # variant_qnames
    sv_info_list = os.popen(''' awk 'NR==%s{print $0}' %s ''' % (svid_x, sv_file_x)).read().strip().split()
    vr_list = sv_info_list[7].split(',')
    with open(vr_txt_x, 'w') as vr_out:
        for vr in vr_list:
            newline = vr + '\n'
            vr_out.write(newline)
    # candidate region bam
    os.system("{samtools} view -@ 10 -h -b {bam_ont} {region} > {bam_out}".format(
        samtools=samtools, bam_ont=ont_bam_x, region=region_x, bam_out=ont_bam_region_x))
    # variant bam
    os.system("({samtools} view -H {inbam}; {samtools} view {inbam}|grep -f {qnames})|"
              "{samtools} view -b - -o {out}".format(
                samtools=samtools, inbam=ont_bam_region_x, qnames=vr_txt_x, out=vr_bam_x))
    os.system('{samtools} index {bam}'.format(samtools=samtools, bam=vr_bam_x))
    # bam to fastq, forward
    os.system(''' %s view -F 2048 %s | awk 'BEGIN {FS="\\t"} {print "@"$1"\\n"$10"\\n+\\n"$11}' - > %s ''' % (
        samtools, vr_bam_x, vr_fq_x))


pools = Pool(5)
for sample_ins in sample_list:
    sample_c = sample_ins + 'C'
    sample_ont = '%s_Cancer_cells' % sample_ins
    svregion_bed = os.path.join(ins_region_dir, "%s_somatic_select_CandidateRegion.bed" % sample_ins)
    if os.path.exists(svregion_bed):
        sample_dir = os.path.join(asmb_dir, sample_c)
        with open(svregion_bed, 'r') as f:
            for line in f:
                c = line.strip().split()
                svid = c[3]
                chrom = c[0]
                svloc = '_'.join([sample_c] + c)
                # qnames
                sv_dir = os.path.join(sample_dir, svloc)
                vr_dir = os.path.join(sv_dir, 'variant_reads')
                if not os.path.exists(vr_dir):
                    os.makedirs(vr_dir)
                sv_file = os.path.join(ins_sniffles_dir, "{sample}/{sample}_INS_sniffles_{chr}_row.bed".format(
                    sample=sample_c, chr=chrom))
                vr_txt = os.path.join(vr_dir, '%s_variant_ont_reads.txt' % svloc)
                # region bam
                region = "%s:%s-%s" % (chrom, int(c[1]) - 500, int(c[1]) + 500)
                ont_bam = "%s/%s/%s.sorted_ngmlr.bam" % (ont_bam_dir, sample_ont, sample_ont)
                ont_bam_region = os.path.join(vr_dir, '%s_candidate_region.bam' % svloc)
                vr_bam = os.path.join(vr_dir, '%s_variant_ont_reads.bam' % svloc)
                vr_fq_x = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
                pools.apply_async(variant_ont_qnames, args=(
                    svid, sv_file, vr_txt, ont_bam, region, ont_bam_region, vr_bam, vr_fq_x))

pools.close()
pools.join()
del pools


########################################################################################
# shasta
########################################################################################
def shasta_run_1(vr_fq_x, shasta_dir_x):
    if os.path.exists(shasta_dir_x):
        os.system("rm -r %s" % shasta_dir_x)
    os.system("{shasta} --input {fa} --Reads.minReadLength 100 --MinHash.minHashIterationCount 100 "
              "--MinHash.minFrequency 2 --Align.minAlignedMarkerCount 100 --MarkerGraph.minCoverage 5 "
              "--assemblyDirectory {outdir}".format(shasta=shasta, fa=vr_fq_x, outdir=shasta_dir_x))


def shasta_run_2(vr_fq_x, shasta_dir_x):
    if os.path.exists(shasta_dir_x):
        os.system("rm -r %s" % shasta_dir_x)
    os.system("{shasta} --input {fa} --Reads.minReadLength 100 --MinHash.minHashIterationCount 100 "
              "--MinHash.minFrequency 2 --Align.minAlignedMarkerCount 100 --MarkerGraph.minCoverage 3 "
              "--assemblyDirectory {outdir}".format(shasta=shasta, fa=vr_fq_x, outdir=shasta_dir_x))


def shasta_run_3(vr_fq_x, shasta_dir_x):
    if os.path.exists(shasta_dir_x):
        os.system("rm -r %s" % shasta_dir_x)
    os.system("{shasta} --input {fa} --Reads.minReadLength 100 --MinHash.minHashIterationCount 100 "
              "--MinHash.minFrequency 2 --Align.minAlignedMarkerCount 100 --MarkerGraph.minCoverage 2 "
              "--assemblyDirectory {outdir}".format(shasta=shasta, fa=vr_fq_x, outdir=shasta_dir_x))


def shasta_run_f(vr_fq_x, shasta_dir_x, sample_c_x, svloc_x):
    shasta_run_1(vr_fq_x, shasta_dir_x)
    shasta_fa = os.path.join(shasta_dir_x, 'Assembly.fasta')
    # shasta low coverage
    if not os.path.getsize(shasta_fa):
        shasta_run_2(vr_fq_x, shasta_dir_x)
        if not os.path.getsize(shasta_fa):
            shasta_run_3(vr_fq_x, shasta_dir_x)
            if not os.path.getsize(shasta_fa):
                return '\t'.join([sample_c_x, svloc_x, 'F']) + '\n'
            else:
                return '\t'.join([sample_c_x, svloc_x, '2']) + '\n'
        else:
            return '\t'.join([sample_c_x, svloc_x, '3']) + '\n'
    else:
        return '\t'.join([sample_c_x, svloc_x, '5']) + '\n'


for sample_ins in sample_list:
    sample_c = sample_ins + 'C'
    svregion_bed = os.path.join(ins_region_dir, "%s_somatic_select_CandidateRegion.bed" % sample_ins)
    if os.path.exists(svregion_bed):
        sample_dir = os.path.join(asmb_dir, sample_c)
        shasta_info_txt = os.path.join(sample_dir, '%s_somatic_INS_shasta_parameter_forward.txt' % sample_c)
        with open(shasta_info_txt, 'w') as out_p:
            with open(svregion_bed, 'r') as f:
                for line in f:
                    c = line.strip().split()
                    svid = c[3]
                    chrom = c[0]
                    svloc = '_'.join([sample_c] + c)
                    sv_dir = os.path.join(sample_dir, svloc)
                    vr_dir = os.path.join(sv_dir, 'variant_reads')
                    vr_fq_x = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
                    shasta_dir = os.path.join(sv_dir, 'shasta')
                    newline = shasta_run_f(vr_fq_x, shasta_dir, sample_c, svloc)
                    out_p.write(newline)
                    out_p.flush()

# contig check
for sample_ins in sample_list:
    sample_c = sample_ins + 'C'
    svregion_bed = os.path.join(ins_region_dir, "%s_somatic_select_CandidateRegion.bed" % sample_ins)
    if os.path.exists(svregion_bed):
        sample_dir = os.path.join(asmb_dir, sample_c)
        shasta_info_txt = os.path.join(sample_dir, '%s_somatic_INS_shasta_contig.txt' % sample_c)
        with open(shasta_info_txt, 'w') as info_out:
            with open(svregion_bed, 'r') as f:
                for line in f:
                    c = line.strip().split()
                    svid = c[3]
                    chrom = c[0]
                    svloc = '_'.join([sample_c] + c)
                    sv_dir = os.path.join(sample_dir, svloc)
                    vr_dir = os.path.join(sv_dir, 'variant_reads')
                    vr_fq_x = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
                    shasta_dir = os.path.join(sv_dir, 'shasta')
                    shasta_fa = os.path.join(shasta_dir, 'Assembly.fasta')
                    if os.path.exists(shasta_fa):
                        if os.path.getsize(shasta_fa):
                            contig1 = int(os.popen("wc -l %s" % shasta_fa).read().strip().split()[0]) / 2
                            if contig1 != 1:
                                info_out.write('\t'.join([svloc, 'run1_Failed:%s' % str(contig1)]) + '\n')
                                shasta_run_2(vr_fq_x, shasta_dir)
                                contig2 = int(os.popen("wc -l %s" % shasta_fa).read().strip().split()[0]) / 2
                                if contig2 != 1:
                                    info_out.write('\t'.join([svloc, 'run2_Failed:%s' % str(contig1)]) + '\n')
                                    shasta_run_3(vr_fq_x, shasta_dir)
                                    contig3 = int(os.popen("wc -l %s" % shasta_fa).read().strip().split()[0]) / 2
                                    if contig3 != 1:
                                        info_out.write('\t'.join([svloc, 'run2_Failed:%s' % str(contig1)]) + '\n')
                        else:
                            info_out.write('\t'.join([svloc, 'empty_Failed']) + '\n')
                    else:
                        info_out.write('\t'.join([svloc, 'not_exist_Failed']) + '\n')

# shasta summary
shasta_summary_txt = os.path.join(asmb_dir, 'xdlab_shasta_summary_check.txt')
summary_out = open(shasta_summary_txt, 'w')
for sample_ins in sample_list:
    sample_c = sample_ins + 'C'
    svregion_bed = os.path.join(ins_region_dir, "%s_somatic_select_CandidateRegion.bed" % sample_ins)
    if os.path.exists(svregion_bed):
        sample_dir = os.path.join(asmb_dir, sample_c)
        with open(svregion_bed, 'r') as f:
            for line in f:
                c = line.strip().split()
                svid = c[3]
                chrom = c[0]
                svloc = '_'.join([sample_c] + c)
                sv_dir = os.path.join(sample_dir, svloc)
                vr_dir = os.path.join(sv_dir, 'variant_reads')
                vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
                shasta_dir = os.path.join(sv_dir, 'shasta')
                shasta_fa = os.path.join(shasta_dir, 'Assembly.fasta')
                if os.path.exists(shasta_fa):
                    if os.path.getsize(shasta_fa):
                        contig1 = int(os.popen("wc -l %s" % shasta_fa).read().strip().split()[0]) / 2
                        if contig1 != 1:
                            summary_out.write('\t'.join([sample_c, svloc, 'Contig', str(contig1)]) + '\n')
                        else:
                            summary_out.write('\t'.join([sample_c, svloc, 'Contig', str(contig1)]) + '\n')
                    else:
                        summary_out.write('\t'.join([sample_c, svloc, 'Empty', '.']) + '\n')
                else:
                    summary_out.write('\t'.join([sample_c, svloc, 'None', '.']) + '\n')
    summary_out.flush()

summary_out.close()


########################################################################################
# racon medaka
########################################################################################
def racon_medaka_run(rc_medaka_dir_x, svloc_x, vr_fq_x, shasta_fa_x, sv_dir_x):
    # run 1
    v2s_sam = os.path.join(rc_medaka_dir_x, '%s_vr_to_shasta.sam' % svloc_x)
    os.system("{minimap2} -ax map-ont -t 10 {ref} {fq} | {samtools} view -hS -q 30 -F 4 > {sam}".format(
        ref=shasta_fa_x, fq=vr_fq_x, sam=v2s_sam, minimap2=minimap2, samtools=samtools))
    racon1_fa = os.path.join(rc_medaka_dir_x, '%s_racon1.fasta' % svloc_x)
    os.system("{racon} -t 1 -w 50  {basecalls} {bam} {draft} > {out}".format(
        racon=racon, basecalls=vr_fq_x, draft=shasta_fa_x, bam=v2s_sam, out=racon1_fa))
    # run 2
    r2r_sam_1 = os.path.join(rc_medaka_dir_x, '%s_vr_to_racon1.sam' % svloc_x)
    os.system("{minimap2} -ax map-ont -t 10 {ref} {fq} | {samtools} view -hS -q 30 -F 4 > {sam}".format(
        ref=racon1_fa, fq=vr_fq_x, sam=r2r_sam_1, minimap2=minimap2, samtools=samtools))
    racon2_fa = os.path.join(rc_medaka_dir_x, '%s_racon2.fasta' % svloc_x)
    os.system("{racon} -t 1 -w 50  {basecalls} {bam} {draft} > {out}".format(
        racon=racon, basecalls=vr_fq_x, draft=racon1_fa, bam=r2r_sam_1, out=racon2_fa))
    # run 3
    r2r_sam_2 = os.path.join(rc_medaka_dir_x, '%s_vr_to_racon2.sam' % svloc_x)
    os.system("{minimap2} -ax map-ont -t 10 {ref} {fq} | {samtools} view -hS -q 30 -F 4 > {sam}".format(
        ref=racon2_fa, fq=vr_fq_x, sam=r2r_sam_2, minimap2=minimap2, samtools=samtools))
    racon3_fa = os.path.join(rc_medaka_dir_x, '%s_racon3.fasta' % svloc_x)
    os.system("{racon} -t 1 -w 50  {basecalls} {bam} {draft} > {out}".format(
        racon=racon, basecalls=vr_fq_x, draft=racon2_fa, bam=r2r_sam_2, out=racon3_fa))
    # run 4
    r2r_sam_3 = os.path.join(rc_medaka_dir_x, '%s_vr_to_racon3.sam' % svloc_x)
    os.system("{minimap2} -ax map-ont -t 10 {ref} {fq} | {samtools} view -hS -q 30 -F 4 > {sam}".format(
        ref=racon3_fa, fq=vr_fq_x, sam=r2r_sam_3, minimap2=minimap2, samtools=samtools))
    racon4_fa = os.path.join(rc_medaka_dir_x, '%s_racon4.fasta' % svloc_x)
    os.system("{racon} -t 1 -w 50  {basecalls} {bam} {draft} > {out}".format(
        racon=racon, basecalls=vr_fq_x, draft=racon3_fa, bam=r2r_sam_3, out=racon4_fa))
    # medaka
    medaka_out = os.path.join(rc_medaka_dir_x, 'racon4_medaka')
    os.system("medaka_consensus -i {basecalls} -d {draft} -o {out} -t 1 -m r941_prom_high_g360".format(
        basecalls=vr_fq_x, draft=racon4_fa, out=medaka_out))
    # to hg38
    medaka_out = os.path.join(rc_medaka_dir_x, 'racon4_medaka')
    rc_medaka_fa = os.path.join(medaka_out, 'consensus.fasta')
    tohg38_dir = os.path.join(sv_dir_x, 'to_hg38')
    if not os.path.exists(tohg38_dir):
        os.makedirs(tohg38_dir)
    rm2r_sam = os.path.join(tohg38_dir, '%s_racon_medaka_to_hg38.sam' % svloc_x)
    os.system("{minimap2} -ax map-ont -t 10 {ref} {fq} > {sam}".format(
        ref=hg38, fq=rc_medaka_fa, sam=rm2r_sam, minimap2=minimap2))
    after_align(rm2r_sam)


pools = Pool(10)
for sample_ins in sample_list[:-1]:
    sample_c = sample_ins + 'C'
    svregion_bed = os.path.join(ins_region_dir, "%s_somatic_select_CandidateRegion.bed" % sample_ins)
    if os.path.exists(svregion_bed):
        sample_dir = os.path.join(asmb_dir, sample_c)
        with open(svregion_bed, 'r') as f:
            for line in f:
                c = line.strip().split()
                svid = c[3]
                chrom = c[0]
                svloc = '_'.join([sample_c] + c)
                sv_dir = os.path.join(sample_dir, svloc)
                vr_dir = os.path.join(sv_dir, 'variant_reads')
                vr_fq = os.path.join(vr_dir, '%s_variant_ont_reads.fastq' % svloc)
                shasta_dir = os.path.join(sv_dir, 'shasta')
                shasta_fa = os.path.join(shasta_dir, 'Assembly.fasta')
                if os.path.exists(shasta_fa):
                    if os.path.getsize(shasta_fa):
                        contig1 = int(os.popen("wc -l %s" % shasta_fa).read().strip().split()[0]) / 2
                        if contig1 == 1:
                            rc_medaka_dir = os.path.join(sv_dir, 'racon_medaka')
                            if os.path.exists(rc_medaka_dir):
                                os.system("rm -r %s" % rc_medaka_dir)
                            os.mkdir(rc_medaka_dir)
                            # racon_medaka_run(rc_medaka_dir, svloc, vr_fq, shasta_fa, sv_dir)
                            pools.apply_async(racon_medaka_run, args=(rc_medaka_dir, svloc, vr_fq, shasta_fa, sv_dir))

pools.close()
pools.join()
del pools
