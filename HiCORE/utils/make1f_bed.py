import sys
import subprocess as sp
import argparse 
import argcomplete

usage = f'\n\t{sys.argv[0]} chrom_size_file restriction_sites_file chromosomes output_dir\n'

parser = argparse.ArgumentParser(description = "", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', dest = 'genomesizeFile', required=True, type=str,  help = 'path to chrom.sizes file')
parser.add_argument('-r', dest = 'restrictionFile', required = True, help = 'path to wholegenome digestion file, space delimited position file')
parser.add_argument('-c', dest = 'chromosomes', default='all', nargs='*', help = 
        '''space-separated chromosome names, or "all".
        Chromosome names must be mathced with the name in the genome size file.''')
parser.add_argument('-o', dest = 'outdir', default = './', help = 'path to output directory')
args = parser.parse_args()
argcomplete.autocomplete(parser)

chrom_size = args.genomesizeFile
res_file = args.restrictionFile
chromosomes = args.chromosomes
out_dir = args.outdir
chrom_file = open(f'{chrom_size}').readlines()
chr_list = []

if chromosomes == 'all':
    for i in chrom_file:
        chr_list.append(i.split('\t')[0].strip())
else:
    chr_list = chromosomes


#a choose certain restriction sites in restriction_sites files
#b column4, numbering 
#c choose the chromosome in restriction_sites files

def makeUnitbed():
    b=1
    c=0
    with open(f'{chrom_size}') as c1, open(f'{res_file}') as r1:
        chrom_context = c1.readlines()
        res_context = r1.readlines()
        writefile = open(f'{out_dir}/1f_unit.bed', 'w')

        for i in chrom_context:
            chr_num = i.split('\t')[0].strip()
            chr_length = i.split('\t')[1].strip()

            frag_sites = res_context[c]
            fragments=frag_sites.split(' ')

            if chr_num in chr_list:
                if chr_num[:-1] != chr_list[0][:-1]:
                    print(f'Chromosome name error### {chr_num[:-1]} != {chr_list[0][:-1]}')
                    sys.exit()
                a=0
                while 1:

                    if a==0:
                        bin_num = str(b)
                        bin_start = str(0)
                        bin_end = str(fragments[a+1].strip())
                        line = chr_num + '\t' + bin_start + '\t' + bin_end + '\t' + bin_num + '\n'
                        writefile.write(line)

                    elif a+1 == len(fragments):
                        break
                    else:
                        bin_num = str(b)
                        bin_start = str(fragments[a].strip())
                        bin_end = str(fragments[a+1].strip())
                        line = chr_num + '\t' + bin_start + '\t' + bin_end + '\t' + bin_num + '\n'
                        writefile.write(line)
                    a += 1
                    b += 1
            c += 1
        writefile.close()

makeUnitbed()
