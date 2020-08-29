import sys
import subprocess as sp
import argparse
import os

#import makeRandomBin

parser = argparse.ArgumentParser(description = 'Extract intra-chromosomal fragment matrices from a ".hic" file and merge them into a single matrix.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', dest = 'hicFile', required=True, help = 'input hic file path')
parser.add_argument('-g', dest = 'genomesizeFile', required=True, help = 'path to chrom.sizes file')
parser.add_argument('-b', dest = 'bedFile', required=True, help = '1f resolution bed file, made by make_1fbed.py')
parser.add_argument('-j', dest = 'juicertools', required=True, help = 'path to juicer_tools.jar file')
parser.add_argument('-d', dest = 'outdir', default = './', help = 'output file directory')
parser.add_argument('-p', dest = 'prefix', default = 'HiCORE_matrix', help = 'output file name prefix, {prefix}.matrix.txt will be created')
parser.add_argument('-c', dest = 'chromosomes', default = "all", nargs='*', help = 'space-delimited specific chromosomes or "all" for all in chrom.sizes file')
args = parser.parse_args()

hic_file = args.hicFile
bed_file = args.bedFile
out_dir = args.outdir
prefix = args.prefix
chromosomes = args.chromosomes
juicer = args.juicertools
chrom_sizes = args.genomesizeFile


## set chromsomes
def makeChromlist():
    global chr_list
    chrom_file = open(f'{chrom_sizes}').readlines()
    chr_list = []
    if chromosomes == 'all':
        for i in chrom_file:
            chr_list.append(i.split('\t')[0].strip())
    else:
        chr_list = chromosomes


def makeOutdir():

    if os.path.isdir(f'{out_dir}') == False:
        os.mkdir(f'{out_dir}')
        os.mkdir(f'{out_dir}/raw')
    else:
        if os.path.isdir(f'{out_dir}/raw') == False:
            os.mkdir(f'{out_dir}/raw')


## convert raw_bin_number to chromosome_corrected_bin_number
def ExtractIF():
    cat_list = []
    for c in chr_list:
        # dump command of juicer_tools
        sp.call(f'java -jar {juicer} dump observed NONE {hic_file} {c} {c} FRAG 1 {out_dir}/raw/{c}_{c}_1f.dump.txt', shell=True, universal_newlines=True)
        # remove the inter-chromosomal interaction & convert the fragment number based on the "unit bed file"
        bin_plus = sp.check_output(f"less {bed_file} | awk '$1==\"{c}\"' | head -1 | cut -f 4", shell=True, universal_newlines=True).strip()
        bin_end = sp.check_output(f"less {bed_file} | awk '$1==\"{c}\"' | tail -1 | cut -f 4", shell=True, universal_newlines=True).strip()
        awk_string = f'print $1+{bin_plus}, $2+{bin_plus}, $3'
        awk_command = f"'{{awk_string}}'"
        sp.call(f"less {out_dir}/raw/{c}_{c}_1f.dump.txt | awk '{{{awk_string}}}' OFS='\t' | awk '$1 <= {bin_end} && $2 <= {bin_end}' OFS='\t'  > {out_dir}/raw/{c}_{c}_1f.dump.converted.txt", shell=True, universal_newlines=True)
        cat_list.append(f'{out_dir}/raw/{c}_{c}_1f.dump.converted.txt')
    cat_join = ' '.join(cat_list)
    # merge each chromosome
    sp.call(f'cat {cat_join} > {out_dir}/{prefix}.matrix.txt', shell=True, universal_newlines=True)

makeChromlist()
makeOutdir()
ExtractIF()
