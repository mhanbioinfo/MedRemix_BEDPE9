# MedRemix_BEDPE9
stand alone MedRemix with BEDPE9 input

## quick start

### install R packages

- install these R packages

```{R}
install.packages("docopt")
install.packages("tidyverse")
install.packages("arrow")
install.packages("S4Vectors")
install.packages("data.table")
install.packages("plyr")
install.packages("doParallel")
install.packages("readr")
install.packages("MASS")
install.packages("flexmix")

BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome")

## manually install 'countreg'
install.packages("distributions3")
install.packages("countreg", repos="http://R-Forge.R-project.org")
```
- download 'BSgenome.Hsapiens.UCSC.hg38' R package and unzip
- https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html
- Source Package = BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz
- note unzipped R package location for next step

### setup run script

- in 'run_MedRemix_bedpe9_H4H.sh'
- edit paths for BEDPE9 input, output directory and BSgenome.Hsapiens.UCSC.hg38 package location

### run toy sample

- 'GSM5067062_2020_6544_human_bedpe.1000lines.bed.gz'
- `bash run_MedRemix_bedpe9_H4H.sh`

## BEDPE9 to bedpe4medremix

- currently filters out chimeric and unaligned reads
- duplicate fragments are kept (i.e. expanded from BEDPE9's read_count column)
- bedpe4medremix's 'mateid', 'paired' and 'mean_mapq' columns have dummy 0 values as these values are not available from BEDPE9, but are not used by medremix anyway

```
## BEDPE9

read1_seqname  read1_start  read1_end  read2_seqname  read2_start  read2_end  read1_strand  read2_strand  read_count
chr2  65125872  65125960  chr2  65125907  65126003  +  -  1
chr18  36253731  36253827  chr18  36253855  36253951  +  -  2
chr11  1446239  1446335  chr11  1446284  1446374  +  -  5
...

## bedpe4medremix

seqname  mateid  paired  start  end  width  mean_mapq
chr2  0  0  65125872  65126003  131  0
chr18  0  0  36253731  36253951  220  0
chr18  0  0  36253731  36253951  220  0
...
```

