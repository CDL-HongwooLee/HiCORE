# HiCORE: Hi-C analysis for identification of core chromatin looping regions with higher resolutions
HiCORE was designed to define high-resolution core interacting fragments with our advanced binning strategy.  


![GitHub_figure](https://user-images.githubusercontent.com/69840555/122363934-513d1880-cf94-11eb-9cce-567fc3a6113d.png)

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

   1. __Multi-layer & multi-fragment binning ---------------------------- HiCORE

   2. Identification of chromatin loops (Fit-HiC2 or other applications) 
   
   3. __Identification of 'CORE' looping regions ---------------------------------- HiCORE

HiCORE pipeline is composed of two python codes, ```hicBinning.py``` and ```HiCORE.py```.  
   
   ```hicBinning.py``` generates multiple layers of multi-fragment bin-arrays and processes the matrix file into Fit-HiC2 input files.   
   ```HiCORE.py``` identifies the interacting regions with maximum detected-layer frequencies from multiple layers of chromatin loop files.
 
 
## Multi-layer & multi-fragment binning - hicBinning.py

Using the genomic structure files (‘chrom.sizes’ and ‘genome digestion’), ```hicBinning.py``` generates individual files for multiple binning layers in .bed format. The fragment matrix is assigned to a file for each binning layer. HiCORE generates ‘fragments’ and ‘interaction’ files, which are necessary for Fit-HiC2 analysis

### Usage

```
usage: hicBinning.py [-h] -i IF_FILE -g GENOMESIZEFILE -r RESFILE -f MIN_LENGTH 
                     [-n BINNING_NUMBER] [-c [CHROMOSOMES [CHROMOSOMES ...]]] 
                     [-j HICORE_DIR] [-o OUTDIR] [-t THREAD] [-k STEP] [-m MEMORY]

Required arguments:

  -i INTERACTION_MATRIX
                        A HiC interaction matrix including single or multiple chromosome(s) interaction data. 
  -g GENOMESIZEFILE     Path to chrom.sizes file
  -r RESTRICTIONFILE    Path to genome digestion file, including chromosome name and restriction position in each raw.
  -f MIN_FRAGLENGTH     The cut-ff bin size. Restriction fragments below the
                        cut-off size will be merged with neighboring fragments
                        
Optional arguments:

  -n BINNING_NUMBER     Number of layers for HiCORE analysis. 
                        (n<=2) : Only forward-reverse binning strategies will be applied. 
                        (n>3) : Add randomly merged bin (n-2) times. 
                        (default : 2)
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        Space-separated chromosome names, or "all". (default : "all")
  -j HICORE_DIR         HiCORE directory path including hicBinning.py, HiCORE.py script files. (default : "./")
  -o OUTDIR             Path to output directory. (default : "./")
  -t THREAD             Number of threads. If t>=n, running time is remarkably reduced but memory-intensive. (default : 1)
  -k STEP               Must be one of ["all","BinningOnly","AfterBinning"]  
                        **BinningOnly : Only For-Rev + Random binning process will be performed, bed files will be provided.  
                        **AfterBinning : Resume the process after random-binning completed. 
                                         All layers in tmp directory will be used for HiCORE analysis. (default : "all")
  -h, --help           show this help message and exit
  
```

### Example run

```python3 hicBinning.py -i ./Interaction_Matrix.txt -g ./reference/hg19.chrom.sizes -r ./reference/hg19_MboI.txt -f 1000 -n 20 -o ./HiCORE_out -t 10 -c Chr1 ```

### Input data
As input files, HiCORE accepts a fragment unit matrix file. For the analysis, HiCORE needs additional files, ‘chrom.sizes’ and ‘genome digestion’ files, which provide information about genome structure. 

##### Matrix file
1.	The fragment unit matrix file is composed of three columns including fragment1 (Column 1), fragment2 (Column 2), and intra-chromosomal interaction frequency (Column 3). Columns 1 and 2 represent the order of (Nth) restriction fragments that are in 3D contact with each other in a given chromosome. For integration of multiple chromosome data, the fragment number of the next chromosome is continued from the last fragment number of the preceding chromosome. The interaction frequency is represented as an integer value of raw interaction counts. You can easily get the matrix file from a '.hic' file using ```utils/dumpMatrix.py``` script.  

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
The ‘chrom.sizes’ file consists of 2 columns (tab-delimited): chromosome names and chromosome sizes.
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
...
```

##### Genome digestion file
3.	The ‘genome digestion’ file includes chromosome name and genomic positions of all restriction fragments in each chromosome in a row. All data in a row are space-delimited, beginning with the chromosome name. The ‘genome digestion’ file can be made using the ‘generate_site_position.py’ Python code in the Juicer package.
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

Each layer of output overlap_bin.bed files is composed of chromosome name, start position, end position and 'overlapped bin-number'.
Unlike the '.bed file', every single fragment is displayed with 'overlapped bin-number'. The fragments, merged into a single bin, have a same bin-number.

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
These files are compatible directly with Fit-HiC2. If you want to use the normalized matrix, you can generate 'bias' vector using 'HiCKRy.py' in Fit-HiC2 packages.

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
Using the multiple layers of multi-fragment resolution chromatin loop files, HiCORE provides chromatin loops which were specified with maximum detected-layer frequencies.

### usage

```
usage: HiCORE.py [-h] -l LOOP_LIST -b BED_LIST -o OVERLAPBED_LIST -p PREFIX
                 [-d OUTPUT_DIR] [-t THREADS] [-x CUTOFF] -u UNIT_BED [-e]

required arguments:

  -h, --help          show this help message and exit
  
  -l LOOP_LIST        loop list file cut by certain criteria, such as q<0.01.
                      Each file must include (chr1, bin1_start, bin2_start, chr2, chr2_start, chr2_end)
                      in their first 6 columns.
                      
  -b BED_LIST         ".bed" list file. The order of list must correct with loop list order.
                      With default setting, bed files are generated in "out_dir/tmp" directory.
                      
  -o OVERLAPBED_LIST  "overlap_bin.bed" list file. The order of list must correct with loop list order.
                      With default setting, "overlap_bin.bed" files are generated in "out_dir/tmp" directory.
                      
  -p PREFIX           Output file name prefix. Output file path = {out_dir}/{prefix}.HiCORE_out.txt
                                          
  -u UNIT_BED         1f_unit.bed file path. The file was generated by hicbinning.py at outdir/tmp directory
                      

optional arguments:

  -d OUTPUT_DIR       Path to output directory. (default = "./")
  
  -t THREADS          Number of threads. 
                      If thread >= layers, the running time is remarkably reduced but memory-intensive. (default = 1)

  -x CUTOFF           The minimum detected-layer frequency of loops. (the number of detected layers / total layers). 0<=x<=1, default = 0.20
  
  -e                  Except the loops with potential noise. 
                      With this option, HiCORE identifies the loops with the local maximum detected-layer frequency more strictly.
                      We strongly recommand to use this option for more stringent analysis.
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

The HiCORE output file contain information about chromatin loops, loop ID and their detected-layer frequency

```
##chr_name1 start1 end1   chr_name2 start2 end2   loopID_1 loopID_2 Frequency
Chr1    3895    4555    Chr1    12486   12934   loop1   local_1 0.82
Chr1    3766    3895    Chr1    6955    7754    loop1   local_2 0.8
Chr1    4555    4645    Chr1    8481    8547    loop1   local_3 0.74
Chr1    4645    4734    Chr1    11956   12192   loop1   local_4 0.54
Chr1    1804    3704    Chr1    5024    5489    loop1   local_5 0.32
Chr1    1804    3545    Chr1    12211   12486   loop2   local_1 0.2
Chr1    4645    4734    Chr1    23699   24149   loop3   local_1 1.0
```
the loopID_1,2 was defined by whether the unprocessed loops(by HiCORE) are connected.
Loops with same loopID_1 and different loopID_2 were connected in unprocessed looping regions, but they would be separated in processing step 
because each of them reveals the maximum detected-layer frequency.


## Utils
### dumpMatrix.py
```dumpMatrix.py``` extract intra-chromosomal fragment matrices from a '.hic' file and merge them into a matrix.
The code uses the 'Dump command' of juicer packages and in-house script for merging each chromosome data.
**the chromosome(s) (-c option) must be equal to the one used in** ```make1f_bed.py```.
```
usage: dumpMatrix.py [-h] -i HICFILE -g GENOMESIZEFILE -b BEDFILE -j
                     JUICERTOOLS [-d OUTDIR] [-p PREFIX]
                     [-c [CHROMOSOMES [CHROMOSOMES ...]]]


  -i HICFILE            Path to a input ".hic" file. 
  -g GENOMESIZEFILE     Path to chr"om.sizes file" 
  -b BEDFILE            1f resolution bed file, made by "make_1fbed.py"
  -j JUICERTOOLS        Path to juicer_tools.jar file
  -d OUTDIR             Output directory (default: ./)
  -p PREFIX             Output file name prefix, {prefix}.matrix.txt will be
                        created (default: HiCORE)
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        Space-delimited specific chromosomes or "all" for all in the chrom.sizes file (default: "all")
```

### make1f_bed.py
```make1f_bed.py``` makes a single-fragment unit bed file which is required for ```dumpMatrix.py``` code.
```
usage: make1f_bed.py [-h] -g GENOMESIZEFILE -r RESTRICTIONFILE
                     [-c [CHROMOSOMES [CHROMOSOMES ...]]] [-o OUTDIR]

  -g GENOMESIZEFILE     Path to chrom.sizes file
  -r RESTRICTIONFILE    Path to wholegenome digestion file, a space-delimited position file
  -c [CHROMOSOMES [CHROMOSOMES ...]]
                        Space-separated chromosome names, or "all" for all in the chrom.sizes file (default: "all")
                        chromosome names must be mathced with the name in the genome size file.
  -o OUTDIR             Path to output directory (default: ./)

```

### makeRandomBin.py
```makeRandomBin.py``` is involoved in ```hicBinning.py``` script for making multi-layer & multi-fragment bed files.
```
usage: makeRandomBin.py [-h] -c GENOMESIZEFILE -f RESTRICTIONFILE -u UNITBED
                        -m MINFRAGLENGTH [-n RANDOMLAYERNUM] [-o OUTDIR]
                        [-t THREAD]

  -g GENOMESIZEFILE   Path to chrom.sizes file 
  -f RESTRICTIONFILE  Path to wholegenome digestion file, a space-delimited position file 
  -u UNITBED          A single fragment resolution bed file. "make1f_bed.py"
                      can make this file. 
  -m MINFRAGLENGTH    The cut-off bin size. Restriction fragments below the
                      cut-off size will be merged with neighboring fragments
  -n RANDOMLAYERNUM   Number of random layers to make. (default: 2)
  -o OUTDIR           Output directory. (default: ./)
  -t THREAD           Number of threads. If t>=n, running time is remarkably
                      reduced but memory-intensive. (default: 1)

```
