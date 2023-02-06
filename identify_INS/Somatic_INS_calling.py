import os
import sys
import pysam

genomelength = {}
ff = open('/nanopore/ref/genome_hg38.bed', 'r')
line = ff.readline()
while line:
    l = line.strip().split('\t')
    genomelength[l[0]] = int(l[2])
    line = ff.readline()
ff.close()
# allreads = tumorbamfile.fetch(reference='%s' % l[0], start=int('%s' % str(int(l[1]-5))),end=int('%s' % str(int(l[2]+5))))
chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
sample = sys.argv[1]
chr=sys.argv[2]

tumorsample = sample + 'T'
bloodsample = sample + 'N'
sv = 'INS'
# chr = 'chr22'

forlengthsnifflesbed='/nanopore/NewSample/'+tumorsample+'/MainChr_'+tumorsample+'_'+sv+'_sniffles.bed'
forlengthtmpbed='/nanopore/NewSample/somatic/INS/'+tumorsample+'/MainChr_'+tumorsample+'_'+sv+'_sniffles.bed'
forlengthoutbed='/nanopore/NewSample/somatic/INS/'+tumorsample+'/Length-MainChr_'+tumorsample+'_INS_sniffles.bed'
length10perrow=int(int(os.popen('wc -l %s'%(forlengthoutbed)).read().strip().split()[0])*0.1)
length10per=os.popen('sort -k7,7n %s|awk \'NR==%s\' '%(forlengthoutbed,str(length10perrow))).read().strip().split('\t')[6]
selectlength10per = divmod(int(length10per), 10)[0] * 10
binnumber = divmod(int(length10per), 10)[0]

tumorbam='/nanopore/NewSample/'+tumorsample+'/'+tumorsample+'.sorted_ngmlr.bam'
bloodbam='/nanopore/NewSample/'+bloodsample+'/'+bloodsample+'.sorted_ngmlr.bam'

tmpbed = '/nanopore/NewSample/somatic/INS/' + tumorsample + '/' + tumorsample + '_' + sv + '_sniffles_'+chr+'.bed'
tmpbedrow='/nanopore/NewSample/somatic/INS/' + tumorsample + '/' + tumorsample + '_' + sv + '_sniffles_'+chr+'_row.bed'
if not os.path.exists(tmpbed):
    os.system('awk -F "\t" \'{if($1=="%s") print $0}\' %s > %s'%(chr,forlengthsnifflesbed,tmpbed))
if not os.path.exists(tmpbedrow):
    os.system('awk \'{print $0"\t"FNR}\' %s > %s'%(tmpbed,tmpbedrow))
tumordir = '/nanopore/NewSample/somatic/INS/' + tumorsample + '/' + chr
blooddir = '/nanopore/NewSample/somatic/INS/' + bloodsample + '/' + chr
if not os.path.exists(tumordir):
    os.mkdir(tumordir)
if not os.path.exists(blooddir):
    os.mkdir(blooddir)
if not os.path.exists(tumordir+'/SpanReads'):
    os.mkdir(tumordir+'/SpanReads')
if not os.path.exists(blooddir+'/SpanReads'):
    os.mkdir(blooddir+'/SpanReads')
if not os.path.exists(tumordir+'/CountBreakpoint'):
    os.mkdir(tumordir+'/CountBreakpoint')
if not os.path.exists(blooddir+'/CountBreakpoint'):
    os.mkdir(blooddir+'/CountBreakpoint')
os.chdir(tumordir)

tumorout = open(tumordir + '/' + sample + 'C_sniffles_' + sv + '_SupportReadsBK_'+chr+'_new.bed', 'w')
bloodout = open(blooddir + '/' + sample + 'B_sniffles_' + sv + '_SupportReadsBK_'+chr+'_new.bed', 'w')
tumorout1 = open(sample + 'C_sniffles_' + sv + '_SupportReadsBK_SupSupport_'+chr+'.bed', 'w')
# tumorout3 = open(tumordir + '/' + sample + 'C_sniffles_' + sv + '_SupportReadsBK_'+chr+'_summary_new.bed', 'w')
tumorout4 = open(tumordir + '/' + sample + 'C_sniffles_' + sv + '_SupportReadsBK_'+chr+'_missINSsite_new.bed', 'w')
i = 1
ff = open(tmpbed, 'r')
line = ff.readline()
while line:
    l = line.strip().split('\t')
    start = l[1]
    end = l[2]
    if l[1] == l[2]:
        end = str(int(l[2]) + 1)
    supportreads = l[7].split(',')
    flag = 0
    selecttumorbam = tumordir + '/SpanReads/' + sample + 'C_' + l[0] + '-' + start + '-' + end + '_sniffles_'+chr+'_' + str(i) + '.bam'
    os.system('/biosoft/samtools-1.9/samtools view -h -b %s %s > %s' % (tumorbam, l[0] + ':' + str(int(start) - 5) + '-' + str(int(end) + 5), selecttumorbam))
    os.system('/biosoft/samtools-1.9/samtools index %s' % (selecttumorbam))
    selectsupporttumorbam=tumordir + '/SpanReads/selectreads-' + sample + 'C_' + l[0] + '-' + start + '-' + end + '_sniffles_'+chr+'_' + str(i) + '.bam'
    Rnames=' -e'.join(supportreads)
    os.system("(/biosoft/samtools-1.9/samtools view -H {bam}; /biosoft/samtools-1.9/samtools view {bam} | grep -e {names}) | /biosoft/samtools-1.9/samtools view -b - -o {out}".format(bam=selecttumorbam, names=Rnames, out=selectsupporttumorbam))
    os.system('/biosoft/samtools-1.9/samtools index %s'%(selectsupporttumorbam))
    selectbloodbam = blooddir + '/SpanReads/' + sample + 'B_' + l[0] + '-' + start + '-' + end + '_sniffles_'+chr+'_' + str(i) + '.bam'
    os.system('/biosoft/samtools-1.9/samtools view -h -b %s %s > %s' % (bloodbam, l[0] + ':' + str(int(start) - 5) + '-' + str(int(end) + 5), selectbloodbam))
    os.system('/biosoft/samtools-1.9/samtools index %s' % (selectbloodbam))
    tumorbamfile = pysam.AlignmentFile(selecttumorbam, 'rb')
    bloodbamfile = pysam.AlignmentFile(selectbloodbam, 'rb')
    INSSTART = []
    ##select range and cut 10bp bin
    tmprangefile = tumordir + '/CountBreakpoint/bin10-' + str(i) + '.bed'
    tmprange = open(tmprangefile, 'w')
    tmpstart = str(int((float(l[1])+float(l[2]))/2))
    tmpend = 0
    num = 1
    for c in range(0, binnumber):
        tmpend = int(tmpstart) + 10
        if tmpend <= genomelength[l[0]]:
            newline = '\t'.join([chr,tmpstart, str(tmpend), str(num)]) + '\n'
            tmpstart = str(tmpend)
            num += 1
            tmprange.write(newline)
            tmprange.flush()
        else:
            newline = '\t'.join([chr,tmpstart, str(genomelength), str(num)]) + '\n'
            tmprange.write(newline)
            tmprange.flush()
            break
    tmpstart = l[1]
    tmpend = 0
    num = 1
    for e in range(0, binnumber):
        tmpend = int(tmpstart) -  10
        if tmpend >= 0:
            newline = '\t'.join([chr,str(tmpend), tmpstart, str(num)]) + '\n'
            tmpstart = str(tmpend)
            num += 1
            tmprange.write(newline)
            tmprange.flush()
        else:
            newline = '\t'.join([chr,'0', tmpstart, str(num)]) + '\n'
            tmprange.write(newline)
            tmprange.flush()
            break
    tmprange.close()
    # scan tumor reads
    infos=tumorbamfile.fetch()
    for info in infos:
        if info.query_name in supportreads:
            # tumorreads.append(sreads)
            if chr_list[info.rname] == l[0] and info.reference_start <= int(start) and info.reference_end >= int(end):
                aa = info.get_reference_positions(full_length=True)
                FindStart = info.qstart
                FindEnd = info.qend
                outstarts = {}
                tmplen = 1
                outstart = 0
                for j in range(FindStart + 1, FindEnd - 1):
                    if str(aa[j]) == str(aa[j - 1]) == 'None':
                        tmplen += 1
                        if str(aa[j + 1]) != 'None' and tmplen >= 50 and str(aa[j + 1]) != 'None':
                            outstarts[str(aa[j + 1] - 1)] = str(tmplen)
                            tmplen = 1
                    else:
                        tmplen = 1
                if outstarts:
                    for z in outstarts.keys():
                        INSSTART.append(int(z))
                else:
                    if int(l[1]) - 10 <= info.reference_end <= int(l[2]) + 10:
                        INSSTART.append(info.reference_end)
                    elif int(l[1]) - 10 <= info.reference_start <= int(l[2]) + 10:
                        INSSTART.append(info.reference_start)
                    else:
                        pass
    ###count each bin support breakpoint
    tmpforRfile = tumordir + '/CountBreakpoint/bin10_count-' + str(i) + '_raw.bed'
    tmpforRfilesort = tumordir + '/CountBreakpoint/bin10_count-' + str(i) + '.bed'
    tmpforR = open(tmpforRfile, 'w')
    if INSSTART:
        instmpbed=tumordir + '/CountBreakpoint/InsPos-' + str(i) + '.bed'
        instmpbedout=open(instmpbed,'w')
        for ins in INSSTART:
            tmpnewline=chr+'\t'+str(ins)+'\t'+str(int(ins)+1)+'\t'+str(i)+'\n'
            instmpbedout.write(tmpnewline)
        instmpbedout.close()
        intersectbed=tumordir + '/CountBreakpoint/bin10-inspos-' + str(i) + '.bed'
        os.system('/software/bedtools2/bedtools2/bin/bedtools intersect -a %s -b %s -c > %s'%(tmprangefile,instmpbed,intersectbed))
        Bin25 = os.popen('awk -F "\t" \'{print $4}\' %s|sort -u' % (tmprangefile)).read().strip().split('\n')
        for Bin in Bin25:
            Bincount=0
            infos = os.popen('awk -F "\t" \'{if($4=="%s") print $5}\' %s' % (Bin, intersectbed)).read().strip().split('\n')
            for info in infos:
                Bincount+=int(info)
            newline='\t'.join([Bin,str(Bincount)])+'\n'
            tmpforR.write(newline)
        os.system('rm %s'%(intersectbed))


    else:
        tumorout4.write(line)
    tmpforR.close()
    os.system('sort -k1,1n %s>%s' % (tmpforRfile, tmpforRfilesort))
    os.system('rm %s' % (tmpforRfile))
    # for R plot to find inflection point
    if os.path.getsize(tmpforRfilesort):
        predictbin = os.popen('Rscript /nanopore/program/select_bin_20200723.R %s' % (tmpforRfilesort)).read().strip().split()[1]
        if predictbin!='NA':
            outrangefile = tumordir + '/SpanReads/' + sample + 'C_' + l[0] + '-' + start + '-' + end + '_selectrange_' + str(i) + '.bed'
            outrange = open(outrangefile, 'w')
            minrange = int((float(l[1])+float(l[2]))/2 - float(predictbin) * 10)
            maxrange = int((float(l[1])+float(l[2]))/2 + float(predictbin) * 10)
            outrangeline='\t'.join([l[0],str(minrange),str(maxrange),str(i)])+'\n'
            outrange.write(outrangeline)
            outrange.close()
            # rescan tumor reads
            flag = 0
            tumorpos = []
            tumorreads = []
            bloodpos = []
            bloodreads = []
            infos = tumorbamfile.fetch()
            for info in infos:
                if info.query_name in supportreads:
                    tumorreads.append(info.query_name)
                    if chr_list[info.rname] == l[0] and info.reference_start <= int(start) and info.reference_end >= int(end):
                        aa = info.get_reference_positions(full_length=True)
                        FindStart = info.qstart
                        FindEnd = info.qend
                        outstarts = {}
                        tmplen = 1
                        outstart = 0
                        for j in range(FindStart + 1, FindEnd - 1):
                            if str(aa[j]) == str(aa[j - 1]) == 'None':
                                tmplen += 1
                                if str(aa[j + 1]) != 'None' and tmplen >= 50 and str(aa[j + 1]) != 'None':
                                    if minrange <= aa[j + 1] - 1 <= maxrange:
                                        outstarts[str(aa[j + 1] - 1)] = str(tmplen)
                                    tmplen = 1
                            else:
                                tmplen = 1
                        if outstarts:
                            for z in outstarts.keys():
                                if int(l[1]) - 25 <= z <= int(l[2]) + 25:
                                    flag = 1
                                    break
                            if flag == 0 and int(l[1]) - 10 <= info.reference_end <= int(l[2]) + 10:
                                outstarts[str(info.reference_end)] = 'SUP'
                            elif flag == 0 and int(l[1]) - 10 <= info.reference_start <= int(l[2]) + 10:
                                outstarts[str(info.reference_start)] = 'SUP'
                            else:
                                pass
                            for k in outstarts.keys():
                                tumorpos.append(int(k))
                                newline = l[0] + '\t' + start + '\t' + end + '\t' + info.query_name + '\t' + k + '\t' + outstarts[k] + '\t' + str(i) + '\t' + l[5] + '\t' + l[6] + '\n'
                                tumorout.write(newline)
                                tumorout.flush()
                        else:
                            if int(l[1]) - 10 <= info.reference_end <= int(l[2]) + 10:
                                tumorpos.append(int(info.reference_end))
                                newline = l[0] + '\t' + start + '\t' + end + '\t' + info.query_name + '\t' + str(info.reference_end) + '\t' + 'SUP' + '\t' + str(i) + '\t' + l[5] + '\t' + l[6] + '\n'
                                tumorout.write(newline)
                                tumorout.flush()
                                # outstarts[str(info.reference_end)]='SUP'
                            elif int(l[1]) - 10 <= info.reference_start <= int(l[2]) + 10:
                                tumorpos.append(int(info.reference_start))
                                newline = l[0] + '\t' + start + '\t' + end + '\t' + info.query_name + '\t' + str(info.reference_start) + '\t' + 'SUP' + '\t' + str(i) + '\t' + l[5] + '\t' + l[6] + '\n'
                                tumorout.write(newline)
                                tumorout.flush()
                                # outstarts[str(info.reference_start)] = 'SUP'
                            else:
                                pass

            # scan blood reads
            allreads = bloodbamfile.fetch(reference='%s' % l[0], start=int('%s' % str(int(start) - 5)),end=int('%s' % str(int(end) + 5)))
            for breads in allreads:
                # if breads.mapping_quality >= 20:
                b_aa = breads.get_reference_positions(full_length=True)
                b_FindStart = breads.qstart
                b_FindEnd = breads.qend
                b_outstarts = {}
                b_tmplen = 1
                b_outstart = 0
                for j in range(b_FindStart + 1, b_FindEnd - 1):
                    if str(b_aa[j]) == str(b_aa[j - 1]) == 'None':
                        b_tmplen += 1
                        if str(b_aa[j + 1]) != 'None' and b_tmplen >= 50 and str(b_aa[j + 1]) != 'None':
                            if minrange <= b_aa[j + 1] - 1 <= maxrange:
                                b_outstarts[str(b_aa[j + 1] - 1)] = str(b_tmplen)
                            b_tmplen = 1
                    else:
                        b_tmplen = 1
                if b_outstarts:
                    bloodreads.append(breads.query_name)
                    for x in b_outstarts.keys():
                        bloodpos.append(int(x))
                        bnewline = l[0] + '\t' + start + '\t' + end + '\t' + breads.query_name + '\t' + x + '\t' + b_outstarts[x] + '\t' + str(i) + '\t' + '.' + '\t' + '.' + '\n'
                        bloodout.write(bnewline)
                        bloodout.flush()
                else:
                    if int(l[1]) - 10 <= breads.reference_end <= int(l[2]) + 10:
                        bloodreads.append(breads.query_name)
                        bloodpos.append(int(breads.reference_end))
                        bnewline = l[0] + '\t' + start + '\t' + end + '\t' + breads.query_name + '\t' + str(breads.reference_end) + '\t' + 'SUP' + '\t' + str(i) + '\t' + '.' + '\t' + '.' + '\n'
                        bloodout.write(bnewline)
                        bloodout.flush()
                    elif int(l[1]) - 10 <= breads.reference_start <= int(l[2]) + 10:
                        bloodreads.append(breads.query_name)
                        bloodpos.append(int(breads.reference_start))
                        bnewline = l[0] + '\t' + start + '\t' + end + '\t' + breads.query_name + '\t' + str(breads.reference_start) + '\t' + 'SUP' + '\t' + str(i) + '\t' + '.' + '\t' + '.' + '\n'
                        bloodout.write(bnewline)
                        bloodout.flush()
                    else:
                        pass

        # if bloodpos:
        #     sumline = '\t'.join([l[0], start, end]) + '\t' + str(min(tumorpos)) + '-' + str(max(tumorpos)) + '\t' + str(len(set(tumorreads))) + '\t' + l[6] + '\t' + str(min(bloodpos)) + '-' + str(max(bloodpos)) + '\t' + str(len(set(bloodreads))) + '\t' + str(i) + '\n'
        #     tumorout3.write(sumline)
        #     tumorout3.flush()
        # else:
        #     sumline = '\t'.join([l[0], start, end]) + '\t' + str(min(tumorpos)) + '-' + str(max(tumorpos)) + '\t' + str(len(set(tumorreads))) + '\t' + l[6] + '\t' + '0' + '-' + '0' + '\t' + '0' + '\t' + str(i) + '\n'
        #     tumorout3.write(sumline)
        #     tumorout3.flush()
    line = ff.readline()
    i = i + 1

ff.close()
tumorout.close()
bloodout.close()
tumorout1.close()
# tumorout2.close()
# tumorout3.close()
tumorout4.close()

######step 2 KS test####

tumor_bkpos=tumordir + '/' + sample + 'C_sniffles_' + sv + '_SupportReadsBK_'+chr+'_new.bed'
blood_bkpos=blooddir + '/' + sample + 'B_sniffles_' + sv + '_SupportReadsBK_'+chr+'_new.bed'

svcount=os.popen('wc -l %s'%(tmpbed)).read().strip().split()[0]
os.system('Rscript /nanopore/program/PosKStest-Lenwilcox_20200725.R {tumorfile} {bloodfile} {svcount} {workdir}'.format(tumorfile=tumor_bkpos,bloodfile=blood_bkpos,svcount=int(svcount),workdir=tumordir))

####step 3 select somatic####

rawbed='/nanopore/NewSample/somatic/INS/' + tumorsample + '/' + tumorsample + '_' + sv + '_sniffles_'+chr+'_row.bed'
KStestresults=tumordir+'/PosKStest-Lenwilcox_pvalue.txt'
bloodemptys=os.popen('awk -F "\t" \'{if($2=="-1" && $3=="-1") print $1}\' %s'%(KStestresults)).read().strip().split('\n')
somaticcleanbed=tumordir+'/'+sample+'C_somatic-BloodClean.bed'
for i in bloodemptys:
    os.system('awk \'$NR=="%s" \' %s >> %s'%(i,rawbed,somaticcleanbed))
bloodnoemptys=os.popen('awk -F "\t" \'{if($2!="-1" && $3!="-1") print $1}\' %s'%(KStestresults)).read().strip().split('\n')
somaticnocleanbed=tumordir+'/'+sample+'C_somatic-BloodNotClean.bed'
for i in bloodnoemptys:
    # os.system('awk \'{if($NR=="%s") print $0}\' %s >> %s'%(i,rawbed,somaticnocleanbed))
    os.system('awk \'NR=="%s"\' %s >> %s' % (i, rawbed,somaticnocleanbed))


