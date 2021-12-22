#!/usr/bin/env python3
import subprocess as sp
import sys
import argparse
import os
import multiprocessing
import time

##usage
parser = argparse.ArgumentParser(description = "", formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-i', dest = 'Interaction_matrix', required = True, help= 'A HiC interaction matrix including single or multiple chromosome(s) interaction data.')
parser.add_argument('-g', dest = 'genomesizeFile', required=True, help = 'Path to chrom.sizes file')
parser.add_argument('-r', dest = 'restrictionFile', required = True, help = 'Path to genome digestion file, including chromosome name and restriction position in each raw')
parser.add_argument('-f', dest = 'minFraglength', required=True, type = int, help = 'The cutoff bin size. Restriction fragments below the cutoff size will be merged with neighboring fragments')
parser.add_argument('-n', dest = 'layerNumber', default = 2, type=int, help = 
'''
Number of layers used for HiCORE analysis.
(n<=2) : Only forward-reverse binning strategies will be applied.  
(n>2) : Add randomly merged bin (n-2) times.
(default = 2)
''')
parser.add_argument('-c', dest = 'chromosomes', default='all', nargs='*', help = 'Space-separated chromosome names, or "all". (default = "all")')
parser.add_argument('-j', dest = 'HiCORE_dir', default = './', help = 'HiCORE directory path including hicBinning.py, HiCORE.py script files. (default = "./")')
parser.add_argument('-o', dest = 'outdir', default = './', help = 'Path to output directory. (default = "./"')
parser.add_argument('-t', dest = 'thread', default = 1, type=int, help = 'Number of threads. If t>=n, the running time is remarkably reduced but memory-intensive. (default = 1)')
parser.add_argument('-k', dest = 'step', default = 'all', type=str, help = 
''' 
Must be one of ["all","BinningOnly","AfterBinning"] 
**BinningOnly : Only For-Rev + Random binning process will be performed, bed files will be provided. 
**AfterBinning : Resume the process after random-binning completed. All layers in tmp directory will be used for HiCORE analysis
(default = "all")
''')
args = parser.parse_args()


## Set arguments & global variables

dump_file = args.Interaction_matrix
out_dir = args.outdir
script_dir = args.HiCORE_dir
chrom_size = args.genomesizeFile
res_file = args.restrictionFile
min_frag = args.minFraglength
n = args.layerNumber
t = args.thread
chromosomes = args.chromosomes
step = args.step

proc_list = list(range(0,t))
chrom_file = open(f'{chrom_size}').readlines()
chr_list = []

if chromosomes == 'all':
    for i in chrom_file:
        chr_list.append(i.split('\t')[0].strip())
else:
    chr_list = chromosomes


##### tmp directory clearing
def tmpClearing():
    if os.path.isdir(f'{out_dir}') == False:
        os.mkdir(f'{out_dir}')
        os.mkdir(f'{out_dir}/tmp')
    elif os.path.isdir(f'{out_dir}/tmp') == True:
        
        if step == "AfterBinning":
            print('\nskip binning process..')
            #global bed_list,overbed_list
            #bed_list = sp.check_output(f'ls {out_dir}/tmp/ | grep bed | grep merged | grep -v overlap', shell=True, universal_newlines=True).strip().split('\n')
            #overbed_list = sp.check_output(f'ls {out_dir}/tmp/ | grep bed | grep merged | grep overlap', shell=True, universal_newlines=True).strip().split('\n')

        elif os.listdir(f'{out_dir}/tmp') != []:
            print(f"\n******{out_dir}/tmp directory is not empty\n")
            sys.exit()

    else:
        os.mkdir(f'{out_dir}/tmp')


##### make single restriction fragment bed files (unit-bed files)
def makeUnitbed():
    b=1
    c=0
    with open(f'{chrom_size}') as c1, open(f'{res_file}') as r1, open(f'{out_dir}/tmp/1f_unit.bed','w') as o1:
        
        res_context = r1.readlines()
        for i in c1:
            chr_num, chr_length = i.strip().split('\t')
                    
            frag_sites = res_context[c]
            fragments=frag_sites.split(' ')
        
            if chr_num in chr_list:
                a=0
                while 1:

                    if a==0:
                        bin_num = str(b)
                        bin_start = str(0)
                        bin_end = str(fragments[a+1].strip())
                        line = chr_num + '\t' + bin_start + '\t' + bin_end + '\t' + bin_num + '\n'
                        o1.write(line)

                    elif a+1 == len(fragments):
                        break
                    else:
                        bin_num = str(b)
                        bin_start = str(fragments[a].strip())
                        bin_end = str(fragments[a+1].strip())
                        line = chr_num + '\t' + bin_start + '\t' + bin_end + '\t' + bin_num + '\n'
                        o1.write(line)

                    a += 1
                    b += 1
            c += 1

##### make merged fragment '.bed' files

def CreateMergedBin():
    if step != "AfterBinning":
        chr_command = ' '.join(chr_list)
        sp.call(f'python3 {script_dir}/utils/makeRandomBin.py -g {chrom_size} -r {res_file} -m {min_frag} -o {out_dir} -n {n} -u {out_dir}/tmp/1f_unit.bed -t {t} -c {chr_command}', shell=True, universal_newlines=True)
    
    ## set list-varibale of bed files.
    global bed_list,overbed_list
    bed_list = sp.check_output(f'ls {out_dir}/tmp/ | grep bed | grep merged | grep -v overlap', shell=True, universal_newlines=True).strip().split('\n')
    overbed_list = sp.check_output(f'ls {out_dir}/tmp/ | grep bed | grep merged | grep overlap', shell=True, universal_newlines=True).strip().split('\n')


##### Conversion of restriction number to exact position of genomes.

def makeBinDic(bed_file, overbed_file):

    front_tag = bed_file.split('/')[-1].split('_')[0]
    overbed_file = f'{front_tag}_merged{min_frag}.overlap_bin.bed'
    
    ## make [original restriction-number : merged restriction-number] & [merged restriction-number : exact position] dictionaries.

    with open(f'{out_dir}/tmp/{bed_file}') as b1, open(f'{out_dir}/tmp/{overbed_file}') as ob1:
        binnum_overbinnum_dic = {}
        overbinnum_mergebin_dic = {}
        c = 1
        
        for o_line in ob1:
            overbinnum = o_line.split('\t')[3].strip()
            binnum_overbinnum_dic[str(c)] = overbinnum
            c += 1
        for b_line in b1:
            binnum = b_line.split('\t')[3].strip()
            mergebin = '\t'.join(b_line.split('\t')[:3])
            overbinnum_mergebin_dic[binnum] = mergebin
    
    ## read a dump file which has consecutive restriction fragment number through the input chromosomes.
            # (The start restriction-number of certain chromosome == the end of just preceding chromosome restriction-number + 1)
    ## Convert the dump file using two dicionaries.

    with open(dump_file) as d1, open(f'{out_dir}/tmp/{front_tag}_merged{min_frag}.tmp','w') as od1:
        for dump_line in d1:
            dump_bin1, dump_bin2, dump_count = dump_line.rstrip().split('\t')
            ## write the converted values
            dump_mergeout_line = overbinnum_mergebin_dic[binnum_overbinnum_dic[dump_bin1]] + '\t' + overbinnum_mergebin_dic[binnum_overbinnum_dic[dump_bin2]] + '\t' + dump_count + '\n'
            
            od1.write(dump_mergeout_line)
        binnum_overbinnum_dic.clear()
        overbinnum_mergebin_dic.clear()

    sp.call(f"less {out_dir}/tmp/{front_tag}_merged{min_frag}.tmp | sort -T {out_dir}/tmp/ -k1,1 -k2,2n -k5,5n | datamash groupby 1,2,3,4,5,6 sum 7 | awk '$1==$4' > {out_dir}/tmp/{front_tag}_merged{min_frag}.txt", shell=True, universal_newlines=True)

##### Make fithic input files

def makeInputFile(input_file):
    
    ## Make '.fragments.txt.gz' files
    sp.call(f"less {out_dir}/tmp/{input_file} | datamash groupby 1,2,3 sum 7 > {out_dir}/tmp/{input_file}_bin1_datamash.tmp", shell=True, universal_newlines=True)
    sp.call(f"less {out_dir}/tmp/{input_file} | sort -T {out_dir}/tmp/ -k4,4 -k5,5n | datamash groupby 4,5,6 sum 7 > {out_dir}/tmp/{input_file}_bin2_datamash.tmp", shell=True, universal_newlines=True)
    sp.call(f"cat {out_dir}/tmp/{input_file}_bin1_datamash.tmp {out_dir}/tmp/{input_file}_bin2_datamash.tmp | sort -T {out_dir}/tmp -k1,1 -k2,2n | datamash groupby 1,2,3 sum 4 | awk '{{print $1, int($2), int(($2+$3)/2), int($4), 1}}' OFS='\t' > {out_dir}/{input_file}.fragments.txt ; gzip {out_dir}/{input_file}.fragments.txt", shell=True, universal_newlines=True)
    
    ## Make '.interaction.txt.gz' files
    awk_command_i = "'{print $1, int(($2+$3)/2), $4, int(($5+$6)/2), $7}' OFS='\t'"
    sp.call(f"less {out_dir}/tmp/{input_file} | awk {awk_command_i} > {out_dir}/{input_file}.interaction.txt ; gzip {out_dir}/{input_file}.interaction.txt", shell=True, universal_newlines=True)

##### whole processes

if __name__ == '__main__':
    tmpClearing()
    makeUnitbed()
    print('create multi-fragment bed files.....')
    CreateMergedBin()
    if step == 'BinningOnly':
        sys.exit()
    print('Assign the interaction matrix to each layer.....')
    if t>1:
        pool = multiprocessing.Pool(processes=t)
        pool.starmap(makeBinDic, zip(bed_list, overbed_list))
        pool.close()
        pool.join()
        input_datalist = sp.check_output(f'ls {out_dir}/tmp | grep {min_frag}.txt', shell=True, universal_newlines=True).strip().split('\n')
        print(f'Make Fit-HiC2 input files.....')
        pool = multiprocessing.Pool(processes=t)
        pool.map(makeInputFile, input_datalist)
        pool.close()
        pool.join()

    else:
        for bed_file,overbed_file in zip(bed_list,overbed_list):
            makeBinDic(bed_file, overbed_file)
        
        input_datalist = sp.check_output(f'ls {out_dir}/tmp | grep {min_frag}.txt', shell=True, universal_newlines=True).strip().split('\n')
        print(f'Make Fit-HiC2 input files.....')

        for input_file in input_datalist:
            makeInputFile(input_file)

    sp.call(f"rm {out_dir}/tmp/*.tmp", shell=True, universal_newlines=True)
    print(f'Fihished.')
