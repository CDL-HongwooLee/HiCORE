# HiCORE
HiCORE was designed to define core interacting regions with higher resolution by the advanced binning strategy.  
Find documentation and examples at ************

![image](https://user-images.githubusercontent.com/69840555/90509331-1f9cca80-e194-11ea-91bb-a3f1b913bfd9.png)

For any question about HiCORE, please contact miso5103@snu.ac.kr.

## Dependencies
python3  
subprocess  
multiprocessing  
Fit-HiC2  
Juicer  

## Installation
Github Install
```
git clone ******
```

## Input data
HiCORE accepts a fragment unit matrix file as an input file. For the analysis, HiCORE needs additional files, ‘chrom.sizes’ and ‘genome digestion’ files, which provide information of genome structure.

#### Matrix file
The fragment unit matrix file is composed of 3 columns including fragment1, fragment2, interaction-frequency information.
The column 1 and 2 represent that the order of (Nth) restriction fragment in given chromosomes. If multiple chromosomes data are used, the fragment number of the next chromosome is continuing from the last fragment number of the preceding chromosome.
```
1       2       38.0
1       3       1.0
2       3       84.0
3       4       8.0
1       5       1.0
2       5       4.0
...     ...     ...
```

#### Chrom.sizes file
The ‘chrom.sizes’ file consists of 2 columns (tab-delimited), including chromosome names and chromosome sizes.
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
...     ...
```

#### Genome digestion file
The ‘genome digestion’ file includes chromosome name and every restriction fragment position of each chromosome in a single raw. Each component of line is space-delimited and the first component indicates a chromosome name. The ‘genome digestion’ file can easily be made by ‘generate_site_position.py’ python code in Juicer packages. (****)
```
Chr1 311 901 1319 1407 1804 3545 ...  
Chr2 1576 1941 2036 2499 2809 3119 ...
...
```

## Run HiCORE
### Whole procedure :

   1. Multi-layers & multi-fragments binning **__HiCORE__

   2. Identification of chromatin loops (Fit-HiC2 or other applications) 
   
   3. Overlap & Expand looping regions **__HiCORE__

HiCORE pipeline is composed of two python codes, hicBinning.py and HiCORE.py. The ```hicBinning.py``` contain codes of multi-fragments binning and processing to Fit-HiC2 input files. The ```HiCORE.py``` contain codes of identifying overlapped & expanded regions of chromatin loop.
 
 
### Multi-layers & multi-fragments binning

Using the genomic structure files (‘chrom.sizes’ file and ‘genome digestion’ file), ```hicBinning.py``` generates multiple layers of multi-fragment binning files in ‘.bed’ formats. Then, the matrix file is assigned to each binning layer and further processed to format of ‘fragments’ and ‘interaction’ files, which are necessary for Fit-HiC2.

#### Usage

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
                        comma-separated chromosome names, or "all". (default : "all")
  -j HICORE_DIR         HiCORE directory path including HiCORE.py script files. (default : "./")
  -o OUTDIR             path to output directory. (default : "./")
  -t THREAD             number of threads, when t>=n, running time is remarkably reduced but memory-intensive. (default : 1)
  -k STEP               Must be one of ["all","BinningOnly","AfterBinning"]  
                        BinningOnly : Only For-Rev + Random binning process will be performed, bed files will be provided.  
                        AfterBinning : Resume the process after random-binning completed, 
                                       all layers in tmp directory will be used for HiCORE analysis. (default : "all")
  -m MEMORY             One of "High" or "Low". "High" means relatively fast process but require more memory.  
                        "Low" process is relatively slow but more stable. (default : "Low")

```

#### Example run

```python3 hicBinning.py -i ./Interaction_Matrix.txt -g ./reference/hg19.chrom.sizes -r ./reference/hg19_DpnII.txt -f 1000 -n 20 -o ./HiCORE_out -t 10 ```

#### Output files

##### .bed file

Each layer of output bed files are composed of chromosome names, start position, end position and bin-numbers.  
Several fragments are merged to single bin. Bed files are generated in the 'out_dir/tmp' directory

```
Chr1    0       1319    1
Chr1    1319    3545    2
Chr1    3545    4555    3
Chr1    4555    5573    4
Chr1    5573    6913    5
```

##### overlap_bin.bed file

Each layer of output overlap_bin.bed files are composed of chromosome names, start position, end position and overlapped bin-numbers.
Unlike '.bed file', every single fragment is displayed with overlapped bin-numbers. The fragments that merged into a single bin have same bin numbers.

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
```
fragments files (gzipped)
```
Chr1    0       659     641     1
Chr1    1319    2432    1304    1
Chr1    3545    4050    3093    1
Chr1    4555    5064    3826    1
Chr1    5573    6243    3263    1
```

### Overlap & Expand looping regions



## Utils
#### mergeDump.py
#### make1f_bed.py
#### makeRandomFragment.py
