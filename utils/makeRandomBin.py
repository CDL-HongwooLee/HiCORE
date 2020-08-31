import sys
import subprocess as sp
import argparse
import random
import multiprocessing
import time

start_time = time.time()

parser = argparse.ArgumentParser(description = 'Make merged FRAG resolution forward/reverse+random-bin abs.bed files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', dest = 'genomesizeFile', required=True, help = 'Path to chrom.sizes file')
parser.add_argument('-r', dest = 'restrictionFile', required=True, help = 'Path to wholegenome digestion file, a space-delimited position file')
parser.add_argument('-u', dest = 'unitBed', required=True, help = 'A single fragment resolution bed file. "make1f_bed.py" can make this file.')
parser.add_argument('-m', dest = 'minFraglength', required=True, type=int, help = 
        '''The cut-off bin size. Restriction fragments below the
        cut-off size will be merged with neighboring fragments''')
parser.add_argument('-n', dest = 'randomLayerNum', default=2, type=int, help = 'Number of random layers to make.')
parser.add_argument('-c', dest = 'chromosomes', default = "all", nargs = '*', help = 'Space-separated chromosome names, or "all"')
parser.add_argument('-o', dest = 'outdir', default = './', help = 'Path to output directory.')
parser.add_argument('-t', dest = 'thread', type=int, default = 1, help = 'Number of threads. If t>=n, running time is remarkably reduced but memory-intensive.')
args = parser.parse_args()

chrom_size = args.genomesizeFile
restriction_sites = open(f'{args.restrictionFile}').readlines()
min_frag = args.minFraglength
out_dir = args.outdir
n = args.randomLayerNum
unit_bed = args.unitBed
t = args.thread
chromosomes = args.chromosomes
chrom_file = open(f'{chrom_size}').readlines()
chr_list = []


## set chromosomes
def makeChromlist():
    global chr_list
    chr_list = []
    if chromosomes == 'all':
        for i in chrom_file:
            chr_list.append(i.split('\t')[0].strip())
    else:
        chr_list = chromosomes

#b column4, numbering 
#c choose the chromosome in restriction_sites files
##forward direction
def Forward_binning():
    """ Merge nearby restriction fragments over specific length from the front end """
    b=1
    c=0
    writefile = open(f'{out_dir}/tmp/Forward_merged{min_frag}.bed','w')
    writefile2 =open(f'{out_dir}/tmp/Forward_merged{min_frag}.overlap_bin.bed','w')

    for i in chrom_file:
        chr_num=i.split('\t')[0].strip()
        chr_length=i.split('\t')[1].strip()
    
        frag_sites=restriction_sites[c]
        fragments=frag_sites.split(' ')

        if chr_num in chr_list:
            a=0
            while a+1<len(fragments):
        
                if a==0:
                    bin_start = 0

                else:
                    bin_start = int(fragments[a].strip())

                bin_end = int(fragments[a+1].strip())
                bin_num = str(b)
                o_bin_line = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                writefile2.write(o_bin_line)

                while 1:

                    if bin_end - bin_start > min_frag:
                        line = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                        writefile.write(line)
                        b=b+1
                        break
                    else:
                        a=a+1
                        try:
                            bin_end = int(fragments[a+1].strip())
                        
                            o_bin_start = int(fragments[a].strip())
                            o_bin_line = chr_num + '\t' + str(o_bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                            writefile2.write(o_bin_line)

                        except IndexError:
                            line = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                            writefile.write(line)
                            b=b+1
                            break
        
                a+=1
        c+=1
    writefile.close()
    writefile2.close()


##reverse direction
def Reverse_binning():
    """ Merge nearby restriction fragments over specific length from the back end to front"""
    
    b=1
    c=0
    writefile3 = open(f'{out_dir}/tmp/Reverse_merged{min_frag}.bed','w')
    writefile4 = open(f'{out_dir}/tmp/Reverse_merged{min_frag}.overlap_bin.bed','w')

    for i in chrom_file:
        chr_num = i.split('\t')[0]
        chr_length = i.split('\t')[1].strip()

        frag_sites = restriction_sites[c]
        fragments=frag_sites.split(' ')
        
        if chr_num in chr_list:

            a=len(fragments)-2

            while a>=0:

                if a==0:
                    bin_start = 0
                    bin_end = fragments[1]
                    bin_num = str(b)
                    line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                    writefile3.write(line2)
                    o_bin_line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                    writefile4.write(o_bin_line2)
                    b=b+1
                    break

                else:
                    bin_start = int(fragments[a].strip())

                bin_end = int(fragments[a+1].strip())
                bin_num = str(b)
                o_bin_line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                writefile4.write(o_bin_line2)
        
                while 1:

                    if bin_end - bin_start > min_frag:
                        line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                        writefile3.write(line2)
                        b=b+1
                        break

                    else:
                        a=a-1
                        try:
                            bin_start = int(fragments[a].strip())

                            o_bin_end = int(fragments[a+1].strip())
                            o_bin_line2 = chr_num + '\t' + str(bin_start) + '\t' + str(o_bin_end) + '\t' + bin_num + '\n'
                            writefile4.write(o_bin_line2)

                        except ValueError:
                            bin_start = 0
                            bin_end = fragments[1]
                            line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                            writefile3.write(line2)
                            o_bin_line2 = chr_num + '\t' + str(bin_start) + '\t' + str(bin_end) + '\t' + bin_num + '\n'
                            writefile4.write(o_bin_line2)

                            b=b+1
                            break
                a=a-1
        c+=1

    writefile3.close()
    writefile4.close()

    sp.call(f'sort -k1,1 -k2,2n -o {out_dir}/tmp/Reverse_merged{min_frag}.bed {out_dir}/tmp/Reverse_merged{min_frag}.bed', shell=True, universal_newlines=True)
    sp.call(f'sort -k1,1 -k2,2n -o {out_dir}/tmp/Reverse_merged{min_frag}.overlap_bin.bed {out_dir}/tmp/Reverse_merged{min_frag}.overlap_bin.bed', shell=True, universal_newlines=True)


## make randomly selected merged bin abs file
random_list = list(range(1,n-1))

def Randombin(trial):
    """ Merge nearby restriction fragments over specific length by random selection strategy """
    random_out_list = []

    with open(f'{unit_bed}') as u1:

        unit_context = u1.readlines()
        res_num = len(unit_context)
        frag_list = list(range(0,res_num))
        binDic = {}
        bin_num = 1

        while 1:

            ran_num = random.choice(frag_list)
            frag_list.remove(ran_num)
            
            frag_bin = unit_context[ran_num].strip()
            chr_name = frag_bin.split('\t')[0]
            frag_length = int(frag_bin.split('\t')[2]) - int(frag_bin.split('\t')[1])
                
            random_out_line = frag_bin + '\t' + str(bin_num) + '\n'
            random_out_list.append(random_out_line)
            binDic[frag_bin] = str(bin_num)
            

            ## add neighbor fragmenst until frag_length > min_frag
            k=1
            m=1
            while frag_length < min_frag:
                if ran_num + k <= res_num - 1 and ran_num -m >= 0:
                    
                    ## add next fragment if the next fragment was not binned before
                    if ran_num + k in frag_list and chr_name == unit_context[ran_num+k].split('\t')[0]:
                        
                        frag_bin2 = unit_context[ran_num + k].strip()
                        frag_length2 = int(frag_bin2.split('\t')[2]) - int(frag_bin2.split('\t')[1])
                        frag_length = frag_length + frag_length2
                        frag_list.remove(ran_num + k)

                        random_out_line = frag_bin2 + '\t' + str(bin_num) + '\n'
                        random_out_list.append(random_out_line)
                        binDic[frag_bin2] = str(bin_num)
                        k = k + 1
                    ## add before fragment if the next fragment was binned
                    elif ran_num - m in frag_list and chr_name == unit_context[ran_num-m].split('\t')[0]:

                        frag_bin2 = unit_context[ran_num - m].strip()
                        frag_length2 = int(frag_bin2.split('\t')[2]) - int(frag_bin2.split('\t')[1])
                        frag_length = frag_length + frag_length2
                        frag_list.remove(ran_num - m)

                        random_out_line = frag_bin2 + '\t' + str(bin_num) + '\n'
                        random_out_list.append(random_out_line)
                        binDic[frag_bin2] = str(bin_num)
                        m = m + 1

                    elif chr_name != unit_context[ran_num + k].split('\t')[0] or chr_name != unit_context[ran_num -m].split('\t')[0]:
                        break

#                    elif frag_list == []:                       
#                        break
                    ## split the jammed multi-fragments and merge to the neighbor fragments
                    else:
                            
                        frag_bin3 = unit_context[ran_num + k].strip()
                        plus_num = binDic[frag_bin3]
                            
                        for i in range(ran_num,ran_num + k):
                            frag_bin2 = unit_context[i].strip()
                            random_out_list.remove(frag_bin2 + '\t' + str(bin_num) + '\n')
                            random_out_list.append(frag_bin2 + '\t' + plus_num + '\n')

                        if m > 1:
                            frag_bin4 = unit_context[ran_num - m].strip()
                            minus_num = binDic[frag_bin4]
                            
                            for j in range(ran_num - m + 1, ran_num):
                                frag_bin2 = unit_context[j].strip()
                                random_out_list.remove(frag_bin2 + '\t' + str(bin_num) + '\n')
                                random_out_list.append(frag_bin2 + '\t' + minus_num + '\n')
                        break

                else:
                    break
            if frag_list == []:
                break
                
            bin_num = bin_num + 1

        out_randcontext = ''.join(random_out_list)
        write_randfile = open(f'{out_dir}/tmp/Random{trial}_merged{min_frag}.tmp','w')
        write_randfile.write(out_randcontext)
        write_randfile.close()
            
        sp.call(f"sort -k5,5n -k4,4n {out_dir}/tmp/Random{trial}_merged{min_frag}.tmp | datamash groupby 1,5 min 2 max 3 | sort -k1,1 -k3,3n | awk '{{print $1, $3, $4, $2}}' OFS='\t' > {out_dir}/tmp/Random{trial}_merged{min_frag}.bed", shell=True, universal_newlines=True)
        sp.call(f"sort -k1,1 -k2,2n {out_dir}/tmp/Random{trial}_merged{min_frag}.tmp | cut -f 1,2,3,5 > {out_dir}/tmp/Random{trial}_merged{min_frag}.overlap_bin.bed", shell=True, universal_newlines=True)

### Process
makeChromlist()
Forward_binning()

if n>=2:
    Reverse_binning()
    if n>2 and t>1:
        pool = multiprocessing.Pool(processes=t)
        pool.map(Randombin,random_list)
        pool.close()

    elif n>2 and t==1:
        for trial in random_list:
            Randombin(trial)
