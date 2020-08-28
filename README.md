# HiCORE
HiCORE was designed to define high-resolution core interacting fragments with our advanced binning strategy.  


![image](https://user-images.githubusercontent.com/69840555/90509331-1f9cca80-e194-11ea-91bb-a3f1b913bfd9.png)


For any question about HiCORE, please contact miso5103@snu.ac.kr.

## Dependencies
python3 >=3.7.  
Fit-HiC2 packages >= 2.0.7 (https://github.com/ay-lab/fithic)
Juicer pacakges (juicer_tools >= 1.13.02) (https://github.com/aidenlab/juicer/wiki)
subprocess, multiprocessing  

## Installation
Github Install
```
git clone https://github.com/CDL-HongwooLee/HiCORE.git
```


# Quick Start
## Whole procedure :

   1. __Multi-layer & multi-fragment binning ---------------------------- HiCORE__

   2. Identification of chromatin loops (Fit-HiC2 or other applications) 
   
   3. __Overlap & Expand looping regions ---------------------------------- HiCORE__

HiCORE pipeline is composed of two python codes, ```hicBinning.py``` and ```HiCORE.py```.  
   
   ```hicBinning.py``` generates multiple layers of multi-fragment bin-arrays and processes the matrix file into Fit-HiC2 input files.   
   ```HiCORE.py``` identifies overlapped & expanded interacting fragments from multiple layers of chromatin loop files.
 
 
## Multi-layer & multi-fragment binning - hicBinning.py

Using the genomic structure files (‘chrom.sizes’ file and ‘genome digestion’ file), ```hicBinning.py``` generates multiple layers of multi-fragment bin-arrays files in ‘.bed’ formats. Then, the matrix file is assigned to each layer and further processed to format of ‘fragments’ and ‘interaction’ files, which are necessary for Fit-HiC2.

### Usage

```
usage: hicBinning.py [-h] -i IF_FILE -g GENOMESIZEFILE -r RESFILE -f MIN_LENGTH 
                     [-n BINNING_NUMBER] [-c [CHROMOSOMES [CHROMOSOMES ...]]] 
                     [-j HICORE_DIR] [-o OUTDIR] [-t THREAD] [-k STEP] [-m MEMORY]

Required arguments:

  -i INTERACTION_MATRIX
                        A HiC interaction matrix including single or multiple chromosome(s) interaction data. 
  -g GENOMESIZEFILE     path to chrom.sizes file
  -r RESTRICTIONFILE    path to genome digestion file, space-delimited position file.
  -f MIN_FRAGLENGTH     The cut-off bin size. restriction fragments below the
                        cut-off size will be merged with neighbouring fragments
                        
Optional arguments:

  -n BINNING_NUMBER     number of layers for HiCORE. if n<=2, only forward-reverse binning strategy will be applied. 
                        n>3, add randomly merged bin (n-2) times. (default : 2)
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        space-separated chromosome names, or "all". (default : "all")
  -j HICORE_DIR         HiCORE directory path including HiCORE.py script files. (default : "./")
  -o OUTDIR             path to output directory. (default : "./")
  -t THREAD             number of threads, when t>=n, running time is remarkably reduced but memory-intensive. (default : 1)
  -k STEP               Must be one of ["all","BinningOnly","AfterBinning"]  
                        BinningOnly : Only For-Rev + Random binning process will be performed, bed files will be provided.  
                        AfterBinning : Resume the process after random-binning completed, 
                                       all layers in tmp directory will be used for HiCORE analysis. (default : "all")
  -m MEMORY             One of "High" or "Low". "High" means relatively fast process but require more memories.  
                        "Low" process is relatively slow but more stable. (default : "Low")
  -h, --help           show this help message and exit
  
```

### Example run

```python3 hicBinning.py -i ./Interaction_Matrix.txt -g ./reference/hg19.chrom.sizes -r ./reference/hg19_DpnII.txt -f 1000 -n 20 -o ./HiCORE_out -t 10 -c Chr1 ```

### Input data
HiCORE accepts a fragment unit matrix file as an input file. For the analysis, HiCORE needs additional files, ‘chrom.sizes’ and ‘genome digestion’ files, which contain information of the genomic structure.

##### Matrix file
The single fragment resolution of matrix file is composed of 3 columns including fragment1, fragment2, interaction-frequency information.
The column 1 and 2 represent that the order of (Nth) restriction fragment in given chromosomes. If multiple-chromosome data are used, the fragment number of the next chromosome is continued from the last fragment number of the preceding chromosome. You can easily get the matrix file from a '.hic' file using ```utils/dumpMatrix.py``` script.  

**Note that chromosomes(-c option) used in all HiCORE scripts including utility scripts, must be identical**

```
1       2       38.0
1       3       1.0
2       3       84.0
3       4       8.0
1       5       1.0
2       5       4.0
...
```

##### Chrom.sizes file
The ‘chrom.sizes’ file consists of 2 columns (tab-delimited), including chromosome names and chromosome sizes.
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
...
```

##### Genome digestion file
The ‘genome digestion’ file includes chromosome name and every restriction fragment position of each chromosome in a single raw. Each component of line is space-delimited and the first component indicates a chromosome name. The ‘genome digestion’ file can easily be made by ‘generate_site_position.py’ python code in Juicer packages. 
```
Chr1 311 901 1319 1407 1804 3545 ...  
Chr2 1576 1941 2036 2499 2809 3119 ...
...
```

### Output files

##### .bed file

Each layer of output bed files is composed of chromosome name, start position, end position and bin-number.  
Several fragments are merged to a single bin. Bed files are generated in the 'out_dir/tmp' directory

```
Chr1    0       1319    1
Chr1    1319    3545    2
Chr1    3545    4555    3
Chr1    4555    5573    4
Chr1    5573    6913    5
...
```

##### overlap_bin.bed file

Each layer of output overlap_bin.bed files is composed of chromosome name, start position, end position and overlapped bin-number.
Unlike '.bed file', every single fragment is displayed with overlapped bin-number. The fragments that merged into a single bin have a same bin-number.

```
Chr1    0       311     1
Chr1    311     901     1
Chr1    901     1319    1
Chr1    1319    1407    2
Chr1    1407    1804    2
Chr1    1804    3545    2
Chr1    3545    3704    3
Chr1    3704    3766    3
Chr1    3766    3895    3
Chr1    3895    4555    3
...
```

##### interaction & fragments file
These files are compatible directly with Fit-HiC2.
If you want to use the normalized matrix, you can generate 'bias' vector using 'HiCKRy.py' in Fit-HiC2 packages.

interaction files (gzipped)
```
Chr1    659     Chr1    659     123
Chr1    659     Chr1    2432    30
Chr1    659     Chr1    4050    13
Chr1    659     Chr1    5064    16
Chr1    659     Chr1    6243    10
...
```
fragments files (gzipped)
```
Chr1    0       659     641     1
Chr1    1319    2432    1304    1
Chr1    3545    4050    3093    1
Chr1    4555    5064    3826    1
Chr1    5573    6243    3263    1
...
```

## Overlap & Expand looping regions -- HiCORE.py
Using the multiple layers of multi-fragment resolution chromatin loop files, HiCORE provides overlapped & expanded interacting fragments.
When N layers are given, HiCORE performes analysis N times, adding layers one by one from the first layer to the Nth layer.

### usage

```
usage: HiCORE.py [-h] -l LOOP_LIST -b BED_LIST -o OVERLAPBED_LIST -p PREFIX
                 [-d OUTPUT_DIR] [-t THREADS]

Required arguments:

  -l LOOP_LIST         loop list file cut by certain criteria, such as q<0.01.
                       Each file must include (chr1, bin1_start, bin2_start, chr2, chr2_start, chr2_end)
                       information in their first 6 columns.
  -b BED_LIST          bed list file. The order of list must correct with loop list order. 
                       At default setting, bed files are generated in "out_dir/tmp" directory.
  -o OVERLAPBED_LIST   overlap_bin.bed list file. The order of list must
                       correct with loop list order. At default setting,
                       overlap_bin.bed files are generated in "out_dir/tmp" directory.
  -p PREFIX            output file name prefix.

Optional arguments:

  -d OUTPUT_DIR        output directory
  -t THREADS           number of threads, when thread >= layers, running time
                       is remarkably reduced. (default = 1)
  -h, --help           show this help message and exit (default = './')
```

### Example run

```python3 /home/HiCORE/HiCORE.py -i ./loop_list.txt -b ./bed_list.txt -o ./overlapbed_list.txt -p Project1 -d ./ -t 10```

### Input data

##### loop_list file
```
Forward/Forward_q0.01.txt
Reverse/Reverse_q0.01.txt
Random1/Random1_q0.01.txt
Random2/Random2_q0.01.txt
Random3/Random3_q0.01.txt
...
```
##### bed_list file
```
tmp/Forward_merged1000.bed
tmp/Reverse_merged1000.bed
tmp/Random1_merged1000.bed
tmp/Random2_merged1000.bed
tmp/Random3_merged1000.bed
...
```
##### overlapbed_list file
```
tmp/Forward_merged1000.overlap_bin.bed
tmp/Reverse_merged1000.overlap_bin.bed
tmp/Random1_merged1000.overlap_bin.bed
tmp/Random2_merged1000.overlap_bin.bed
tmp/Random3_merged1000.overlap_bin.bed
...
```

### Output data
With N number of output files, HiCORE presents recommendable layer-number and their average bin sizes.

   First minima : The layer-number when the bin sizes have first minimum values.  
   Minimum binSize : The layer-number which hava minimum average bin sizes.  
   Whole layers : All given layers are used.  
   
The final output files consist of overlapped/expanded interacting regions (first 6 columns) and original interacting regions (last 6 columns, the loops at first layer).


```
Chr1    3766    3895    Chr1    20066   21081   Chr1    3545    4555    Chr1    19737   21081
Chr1    3766    3895    Chr1    21081   21206   Chr1    3545    4555    Chr1    21081   22236
Chr1    3766    3895    Chr1    23699   24691   Chr1    3545    4555    Chr1    23597   24691
Chr1    3766    3895    Chr1    27471   28274   Chr1    3545    4555    Chr1    27471   29385
Chr1    3766    4555    Chr1    7508    8037    Chr1    3545    4555    Chr1    6913    8037
```

## Utils
### dumpMatrix.py
```dumpMatrix.py``` extract intra-chromosomal fragment matrices from a '.hic' file and merge them into a matrix.
The code uses the 'Dump command' of juicer packages and in-house script for merging each chromosome data.
**the chromosome(s) (-c option) must be equal to the one used in** ```make1f_bed.py```.
```
usage: dumpMatrix.py [-h] -i HICFILE -g GENOMESIZEFILE -b BEDFILE -j JUICERDIR
                     [-o OUTDIR] [-c [CHROMOSOMES [CHROMOSOMES ...]]]

  -i HICFILE            input hic file path 
  -g GENOMESIZEFILE     path to chrom.sizes file 
  -b BEDFILE            1f resolution bed file, made by make_1fbed.py
  -j JUICERTOOLS        path to juicer_tools.jar file
  -o OUTDIR             dump output file directory (default: ./)
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        space-delimited specific chromosomes or "all" in chrom.sizes file (default: "all")
```

### make1f_bed.py
```make1f_bed.py``` makes a single-fragment unit bed file which is required for ```dumpMatrix.py``` code.
```
usage: make1f_bed.py [-h] -g GENOMESIZEFILE -r RESTRICTIONFILE
                     [-c [CHROMOSOMES [CHROMOSOMES ...]]] [-o OUTDIR]

  -g GENOMESIZEFILE     path to chrom.sizes file
  -r RESTRICTIONFILE    path to wholegenome digestion file, space delimited
                        position file
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        space-separated chromosome names, or "all" (default: "all")
                        chromosome names must be mathced with the name in the genome size file.
  -o OUTDIR             path to output directory (default: ./)

```

### makeRandomBin.py
```makeRandomBin.py``` is involoved in ```hicBinning.py``` script for making multi-layer & multi-fragment bed files.
```
usage: makeRandomBin.py [-h] -c GENOMESIZEFILE -f RESTRICTIONFILE -u UNITBED
                        -m MINFRAGLENGTH [-n RANDOMLAYERNUM] [-o OUTDIR]
                        [-t THREAD]

  -g GENOMESIZEFILE   path to chrom.sizes file 
  -f RESTRICTIONFILE  path to chrom.sizes file 
  -u UNITBED          a single fragment resolution bed file. "make1f_bed.py"
                      can make this file. 
  -m MINFRAGLENGTH    The cut-off bin size. restriction fragments below the
                      cut-off size will be merged with neighbouring fragments
  -n RANDOMLAYERNUM   Number of random layers to make. (default: 2)
  -o OUTDIR           path to output directory. (default: ./)
  -t THREAD           number of threads, when t>=n, running time is remarkably
                      reduced but memory-intensive. (default: 1)

```
