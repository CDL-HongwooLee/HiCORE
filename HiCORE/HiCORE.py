import sys
import subprocess as sp
import argparse
import os
import time
import multiprocessing


parser = argparse.ArgumentParser(description = '', formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument('-l', dest = 'loop_list', required=True, type = str, help = 
'''loop list file cut by certain criteria, such as q<0.01.
Each file must include (chr1, bin1_start, bin2_start, chr2, chr2_start, chr2_end)
in their first 6 columns.

''')
parser.add_argument('-b', dest = 'bed_list', required=True, type = str, help = 
'''".bed" list file. The order of list must correct with loop list order.
At default setting, bed files are generated in "out_dir/tmp" directory.

''')
parser.add_argument('-o', dest = 'overlapbed_list', required=True, type = str, help = 
'''"overlap_bin.bed" list file. The order of list must correct with loop list order.
At default setting, "overlap_bin.bed" files are generated in "out_dir/tmp" directory.

''')
parser.add_argument('-p', dest = 'prefix', required=True, help = 'Output file name prefix. Output file path = {out_dir}/{prefix}-{overlap/expand}.{layers}.HiCORE.txt')
parser.add_argument('-d', dest = 'output_dir', default = './', help = 'Path to output directory. (default = "./")')
parser.add_argument('-t', dest = 'threads', default = 1, type=int, help = 
'''Number of threads. 
If thread >= layers, the running time is remarkably reduced but memory-intensive. (default = 1)
''')
args = parser.parse_args()


## Variable setting
bed_list_file = args.bed_list
loop_list_file = args.loop_list
overbed_list_file = args.overlapbed_list
out_dir = args.output_dir
out_file = args.prefix
t = args.threads

split_loop_list = []
binSizelist = []

with open(f'{bed_list_file}') as b1:
    bed_list = b1.readlines()
    first_layer_bed = bed_list[0].strip()
with open(f'{loop_list_file}') as l1:
    loop_list = l1.readlines()
    layer_num = len(loop_list)
    split_loop_list = []
    for filename in loop_list:
        split_loop_file = filename.strip() + '.splitbin.txt'
        split_loop_list.append(split_loop_file)
with open(f'{overbed_list_file}') as o1:
    overbed_list = o1.readlines()
    first_layer_overbed = overbed_list[0].strip()


def tmpClearing():
    if os.path.isdir(f'{out_dir}') == True:
        if os.path.isdir(f'{out_dir}/tmp') == False:
            os.mkdir(f'{out_dir}/tmp')
        elif os.listdir(f'{out_dir}/tmp') != []:
            print(f"\n******{out_dir}/tmp directory is not empty\n")
            sys.exit()
    else:
        os.mkdir(f'{out_dir}')
        os.mkdir(f'{out_dir}/tmp')


## Split the multi-fragments resolution loop files into single fragment resolution
def splitloop(bed_file, overbed_file, loop_file):
    mid_to_binnum = {}
    front_tag = '.'.join(loop_file.split('.')[:-1])

    ### mid point of loop : bin_num dictionary creation
    with open(f'{bed_file}'.strip()) as b1:
        bed_context = b1.readlines()

        for blines in bed_context:
    
            bsplit = blines.split('\t')
            chr_name = bsplit[0]
            mid_point = str(int((int(bsplit[1])+int(bsplit[2]))/2))
            bin_num = bsplit[3].strip()
            mid_to_binnum[chr_name + '\t' + mid_point] = bin_num

    ### (a multi-fragments bin number) : (intact bin of split fragments) dictionary  
    with open(f'{overbed_file}'.strip()) as o1:
        ocontext = o1.readlines()
        binnum_to_fulldic = {}
        obin_list = []

        a=0
        while a <= len(ocontext)-1:
            osplit = ocontext[a].split('\t')
            ofull_bin = '\t'.join(ocontext[a].split('\t')[:3])

            if a == 0:
                obin_list=[]
                obin_list.append(ofull_bin)
                obin_num = str(1)
            else:
                if osplit[3].strip() == obin_num:
                    obin_list.append(ofull_bin)
                    if a == len(ocontext)-1:
                        binnum_to_fulldic[obin_num] = obin_list
                        break
                else:
                    binnum_to_fulldic[obin_num] = obin_list
                    obin_list=[]
                    obin_num = osplit[3].strip()
                    obin_list.append(ofull_bin)
            
            a += 1

    with open(f'{loop_file}'.strip()) as l1:
        loop_context = l1.readlines()
        filename = loop_file.strip()
        writefinal = open(f'{filename}.splitbin.txt','w')

        for llines in loop_context:
            bin1 = '\t'.join(llines.split('\t')[:2])
            bin2 = '\t'.join(llines.split('\t')[2:4])

            full_bin1 = mid_to_binnum[bin1]
            full_bin2 = mid_to_binnum[bin2]

            for i in binnum_to_fulldic[full_bin1]:
                for j in binnum_to_fulldic[full_bin2]:
                    full_bin_out = i + '\t' + j + '\n'
                    writefinal.write(full_bin_out)

        writefinal.close()

## detect expanded & overlapped fragments
def expandloop():
    cat_loop = ' '.join(split_loop_list)
    sp.call(f'cat {cat_loop} | sort -k1,1 -k2,2n -k5,5n | uniq > {out_dir}/tmp/{out_file}-expand.tmp', shell=True, universal_newlines=True)

def detectOverlap(k):
    first_layer = split_loop_list[0]
    if k==1:
        sp.call(f"cp {first_layer} {out_dir}/tmp/{out_file}-overlap.n{k}.tmp", shell=True, universal_newlines=True)
    else:
        grep_command = '| grep -f '.join(split_loop_list[:k])
        k_layer = split_loop_list[k-1]
        sp.call(f"less {out_dir}/tmp/{out_file}-overlap.n{k-1}.tmp | grep -f {k_layer} > {out_dir}/tmp/{out_file}-overlap.n{k}.tmp", shell=True, universal_newlines=True)

## re-merge process
def addBin():
    ### Dictionary for re-merging process using the first layer.
    global bin_to_binnum_accDic
    global binnum_acc_to_binDic
    global bedDic
    with open(f'{first_layer_overbed}') as ob1:
        overbed_context = ob1.readlines()
        overbed_name = '.'.join(first_layer_overbed.split('/')[-1].split('.')[:-1])
        
        overbed_acc_outfile = open(f'{out_dir}/tmp/{out_file}-{overbed_name}.addBin.bed','w')
        overbed_binnum = 0
        bin_acc = 0
        bin_to_binnum_accDic = {}
        for overbed_lines in overbed_context:
            if overbed_binnum == overbed_lines.split('\t')[3].strip():
                bin_acc += 1
            else:
                bin_acc = 1
            overbed_binnum = overbed_lines.split('\t')[3].strip()
            out_lines = overbed_lines.strip() + '\t' + str(bin_acc) + '\n'
            overbed_acc_outfile.write(out_lines)

            bintext = '\t'.join(overbed_lines.split('\t')[:3])
            binnum_acc = overbed_binnum + '\t' + str(bin_acc)
            bin_to_binnum_accDic[bintext] = binnum_acc
            
        binnum_acc_to_binDic = dict((value,key) for key,value in bin_to_binnum_accDic.items())
        overbed_acc_outfile.close()

    with open(f'{first_layer_bed}') as b1:
        bed_context = b1.readlines()
        bedDic = {}

        for bed_lines in bed_context:
            bed_bin = '\t'.join(bed_lines.split('\t')[:3])
            bed_num = bed_lines.split('\t')[3].strip()
            bedDic[bed_num] = bed_bin

def addBinLoop(prefix):
    ### merge the expanded & overlapped loop files into multi-fragment resolution using some dictionaries.
    with open(f'{prefix}.tmp') as a1:
        overlap_context = a1.readlines()
        overlap_outfile = open(f'{prefix}_addBin.tmp','w')      

        for over_line in overlap_context:
            overbin1 = '\t'.join(over_line.split('\t')[:3]).strip()
            overbin2 = '\t'.join(over_line.split('\t')[3:]).strip()

            binacc1 = bin_to_binnum_accDic[overbin1]
            binacc2 = bin_to_binnum_accDic[overbin2]

            overlap_outline = binacc1 + '\t' + binacc2 + '\n'
            overlap_outfile.write(overlap_outline)
        
        overlap_outfile.close()
        sp.call(f'less {prefix}_addBin.tmp | sort -k1,1n -k3,3n -k2,2n -k4,4n -T {out_dir}/tmp/ > {prefix}_sorted.addBin.tmp', shell=True, universal_newlines=True)


def bin2_merge(prefix):
    x=1
    with open(f'{prefix}_sorted.addBin.tmp') as f1:
        context1 = f1.readlines()
        bin2_merged_writefile = open(f'{prefix}_bin2_merged.tmp','w')

        while x <= len(context1)-1:
            split_n0 = context1[x-1].split('\t')
            split_n1 = context1[x].split('\t')

            bin1_ij_n0 = '\t'.join(split_n0[:2])
            bin2_i_n0 = split_n0[2]
            bin2_j_n0 = split_n0[3].strip()
            bin1_ij_n1 = '\t'.join(split_n1[:2])
            bin2_i_n1 = split_n1[2]
            bin2_j_n1 = split_n1[3].strip()

            if bin1_ij_n0 == bin1_ij_n1 and bin2_i_n0 == bin2_i_n1:

                bin2_chr = binnum_acc_to_binDic[bin2_i_n0 + '\t' + bin2_j_n0].split('\t')[0]
                bin2_start = binnum_acc_to_binDic[bin2_i_n0 + '\t' + bin2_j_n0].split('\t')[1]

                while 1:
                    if x == len(context1):
                        bin2_end = binnum_acc_to_binDic[bin2_i_n1 + '\t' + bin2_j_n1].split('\t')[2]
                        bin2_merged_outline = bin1_ij_n0 + '\t' + bin2_chr + '\t' + bin2_start + '\t' + bin2_end + '\t' + bin2_i_n1 + '\n'
                        bin2_merged_writefile.write(bin2_merged_outline)
                        break

                    split_n0 = context1[x-1].split('\t')
                    split_n1 = context1[x].split('\t')
                    bin1_ij_n0 = '\t'.join(split_n0[:2])
                    bin2_i_n0 = split_n0[2]
                    bin2_j_n0 = split_n0[3].strip()
                    bin2_i_n1 = split_n1[2]
                    bin2_j_n1 = split_n1[3].strip()

                    if int(bin2_j_n1) - int(bin2_j_n0) == 1 and bin2_i_n0==bin2_i_n1:
                        x=x+1

                    else:
                        bin2_end = binnum_acc_to_binDic[bin2_i_n0 + '\t' + bin2_j_n0].split('\t')[2]
                        bin2_merged_outline = bin1_ij_n0 + '\t' + bin2_chr + '\t' + bin2_start + '\t' + bin2_end + '\t' + bin2_i_n0 + '\n'
                        bin2_merged_writefile.write(bin2_merged_outline)
                        break
            else:

                bin2_merged_outline = bin1_ij_n0 + '\t' + binnum_acc_to_binDic[bin2_i_n0 + '\t' + bin2_j_n0] + '\t' + bin2_i_n0 + '\n'
                bin2_merged_writefile.write(bin2_merged_outline)

            x+=1

        bin2_merged_writefile.close()
        sp.call(f'less {prefix}_bin2_merged.tmp | sort -k3,3 -k4,4n -k1,1n -k2,2n -T {out_dir}/tmp/> {prefix}_bin2_sorted.tmp', shell=True, universal_newlines=True)


def bin1_merge(prefix):
    y=1

    with open(f'{prefix}_bin2_sorted.tmp') as f2:
        bin1_merged_writefile = open(f'{prefix}_merged_output.tmp', 'w')
        context2 = f2.readlines()

        while y <= len(context2)-1:
            split_n0 = context2[y-1].split('\t')
            split_n1 = context2[y].split('\t')

            bin1_i_n0 = split_n0[0]
            bin1_j_n0 = split_n0[1]
            bin2_ij_n0 = '\t'.join(split_n0[2:5])
            bin1_i_n1 = split_n1[0]
            bin1_j_n1 = split_n1[1]
            bin2_ij_n1 = '\t'.join(split_n1[2:5])

            bin1_binnum = split_n0[0]
            bin2_binnum = split_n0[5].strip()

            bin1_ori = bedDic[bin1_binnum]
            bin2_ori = bedDic[bin2_binnum]

            if bin2_ij_n0 == bin2_ij_n1 and bin1_i_n0 == bin1_i_n1:
                bin1_chr = binnum_acc_to_binDic[bin1_i_n0 + '\t' + bin1_j_n0].split('\t')[0]
                bin1_start = binnum_acc_to_binDic[bin1_i_n0 + '\t' + bin1_j_n0].split('\t')[1]

                while 1:
                    if y == len(context2):
                        bin1_end = binnum_acc_to_binDic[bin1_i_n1 + '\t' + bin1_j_n1].split('\t')[2]
                        bin1_merged_outline = bin1_chr + '\t' + bin1_start + '\t' + bin1_end + '\t' + bin2_ij_n1 + '\t' + bin1_ori + '\t' + bin2_ori + '\n'
                        bin1_merged_writefile.write(bin1_merged_outline)
                        break

                    split_n0 = context2[y-1].split('\t')
                    split_n1 = context2[y].split('\t')

                    bin1_i_n0 = split_n0[0]
                    bin1_j_n0 = split_n0[1]
                    bin2_ij_n0 = '\t'.join(split_n0[2:5])
                    bin1_i_n1 = split_n1[0]
                    bin1_j_n1 = split_n1[1]
                    bin2_ij_n1 = '\t'.join(split_n1[2:5])
                    bin2_binnum2 = split_n1[5].strip()

                    if int(bin1_j_n1) - int(bin1_j_n0) == 1 and bin1_i_n0 == bin1_i_n1 and bin2_binnum == bin2_binnum2:
                        y=y+1

                    else:
                        bin1_end = binnum_acc_to_binDic[bin1_i_n0 + '\t' + bin1_j_n0].split('\t')[2]
                        bin1_merged_outline = bin1_chr + '\t' + bin1_start + '\t' + bin1_end + '\t' + bin2_ij_n0 + '\t' + bin1_ori + '\t' + bin2_ori + '\n'
                        bin1_merged_writefile.write(bin1_merged_outline)
                        break

            else:
                bin1_merged_outline = binnum_acc_to_binDic[bin1_i_n0 + '\t' + bin1_j_n0] + '\t' + bin2_ij_n0 + '\t' + bin1_ori + '\t' + bin2_ori + '\n'
                bin1_merged_writefile.write(bin1_merged_outline)

            y+=1
        
        bin1_merged_writefile.close()
    sp.call(f'less {prefix}_merged_output.tmp | sort -k1,1 -k4,4 -k2,2n -k3,3n -k5,5n -k6,6n -T {out_dir}/tmp/ > {prefix}.HiCORE.txt', shell=True, universal_newlines=True)

def reMergeBin(prefix):
    addBinLoop(prefix)
    bin2_merge(prefix)
    bin1_merge(prefix)


## Post-merge process. 
def evaluateBinsize(k):
    ### Calculate the average bin sizes of each output.
    original_loop_number = sp.check_output(f'less {out_dir}/tmp/{out_file}-overlap.n{k}.HiCORE.txt | cut -f 7- | sort -k1,1 -k2,2n -k5,5n | uniq | wc -l', shell=True, universal_newlines=True).strip()
    sumBinsizes = sp.check_output(f"less {out_dir}/tmp/{out_file}-overlap.n{k}.HiCORE.txt | awk '{{sum += $3-$2 + $6-$5}} END {{print int(sum/2)}}'", shell=True, universal_newlines=True).strip()
    avgBinsizes = int(float(sumBinsizes)/float(original_loop_number))
    binSizelist.append(avgBinsizes)

def printEvaluation():
    ## Print out the recommendable layer numbers.
    i=0
    while 1:
        if i==len(binSizelist) -1:
            First_minima = 'None'
            First_minima_layer = 'None'
            break
        elif binSizelist[i] < binSizelist[i+1]:
            First_minima = str(binSizelist[i+1]) + 'bp'
            First_minima_layer = i + 1
            break
        i=i+1
    
    minimum = min(binSizelist)
    minimum_layer = int(binSizelist.index(minimum)) + 1
    whole_layers = binSizelist[-1]
    print('\n')
    print(f'\tFirst minima :  \tlayer N = {First_minima_layer}, Average bin sizes = {First_minima}')
    print(f'\tMinimum binSize :\tlayer N = {minimum_layer}, Average bin sizes = {minimum}bp')
    print(f'\tWhole layers :  \tlayer N = {layer_num}, Average bin sizes = {whole_layers}bp')
    print('\n')
    sp.call(f'mv {out_dir}/tmp/{out_file}-expand.HiCORE.txt {out_dir}', shell=True, universal_newlines=True)
    if First_minima != 'None':
        sp.call(f'mv {out_dir}/tmp/{out_file}-overlap.n{First_minima_layer}.HiCORE.txt {out_dir}', shell=True, universal_newlines=True)
    if minimum_layer != First_minima:
        sp.call(f'mv {out_dir}/tmp/{out_file}-overlap.n{minimum_layer}.HiCORE.txt {out_dir}', shell=True, universal_newlines=True)
    if layer_num != minimum_layer:
        sp.call(f'mv {out_dir}/tmp/{out_file}-overlap.n{layer_num}.HiCORE.txt {out_dir}', shell=True, universal_newlines=True)

    sp.call(f'mv {out_dir}/tmp {out_dir}/Other_layer', shell=True, universal_newlines=True)

if __name__ == '__main__':
    tmpClearing()
    
    if t>1:
        pool = multiprocessing.Pool(processes=t)
        pool.starmap(splitloop, zip(bed_list, overbed_list, loop_list))
        pool.close()
        pool.join()
    else:
        for bed_file, overbed_file, loop_file in zip(bed_list, overbed_list, loop_list):
            splitloop(bed_file, overbed_file, loop_file)
    
    
    expandloop()
    klist = list(range(1,layer_num + 1))
    for k in klist:
        detectOverlap(k)

    addBin()
    tmp_list = sp.check_output(f'ls {out_dir}/tmp/*.tmp', shell=True, universal_newlines=True).strip().split('\n')
    prefix_list = []
    for tmp_file in tmp_list:
        prefix_list.append('.'.join(tmp_file.split('.')[:-1]))
    
    if t>1:
        pool = multiprocessing.Pool(processes=t)
        pool.map(reMergeBin, prefix_list)
        pool.close()
        pool.join()
    else:
        for prefix in prefix_list:
            reMergeBin(prefix)
    sp.call(f'rm {out_dir}/tmp/*.tmp', shell=True, universal_newlines=True)
    

    for k in klist:
        evaluateBinsize(k)
    printEvaluation()
