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
...   ...
```

#### Genome digestion file
The ‘genome digestion’ file includes chromosome name and every restriction fragment position of each chromosome in a single raw. Each component of line is space-delimited and the first component indicates a chromosome name. The ‘genome digestion’ file can easily be made by ‘generate_site_position.py’ python code in Juicer packages. (****)
```
Chr1 311 901 1319 1407 1804 3545 ...  
Chr2 1576 1941 2036 2499 2809 3119 ...
...
```

## Run HiCORE
#### hicBinning.py
#### HiCORE.py

## Output data

## Utils
#### mergeDump.py
#### make1f_bed.py
#### makeRandomFragment.py
