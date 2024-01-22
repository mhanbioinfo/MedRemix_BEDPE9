suppressMessages(library(docopt))

' bedpe9_for_medremixBEDPE v1.0
Convert BEDPE9 format to input accepted by medremixBEDPE

Usage:
    bedpe9_for_medremixBEDPE.R -i INPUT -o OUT_FILE

Options:
    -i --input INPUT        Path to cleaned BEDPE4MEDREMIX file (full path)
    -o --outfile OUT_FILE   Path to output bin_stat file (full path)
' -> doc

args <- docopt(doc, version='1.0')

suppressMessages(library(tidyverse))
suppressMessages(library(S4Vectors))
suppressMessages(library(data.table))

# args = list()
# bedpe9_fpath = "/Users/minghan/bioinfoproj/PughLabPMH/_projects/MedRemix_bedpe9/GSM5067062_2020_6544_human_bedpe.1000lines.bed.gz"
# bedpe9_fpath = "/Users/minghan/bioinfoproj/PughLabPMH/_projects/MedRemix_bedpe9/GSM5067062_2020_6544_human_bedpe.bed.gz"
# output_fpath = "/Users/minghan/bioinfoproj/PughLabPMH/_projects/MedRemix_bedpe9/outputs/GSM5067062_2020_6544_human_bedpe.1000lines.bedpe4medremix"

bedpe9_fpath = args[['input']]
output_fpath = args[['outfile']]

message('Reading in bedpe.gz file...')
bedpe9_f = data.table::fread(bedpe9_fpath)
# bedpe9_f %>% head()
# bedpe9_f %>% names()

bedpe9_f.df = DataFrame(bedpe9_f)
bedpe9_f.df.dup_expanded = rep(bedpe9_f.df, bedpe9_f.df$V9)
# bedpe9_f.df.dup_expanded %>% head()

message('Converting BEDPE9 into format accepted by MedRemixBEDPE...')
bedpe9_for_medremix =
  bedpe9_f.df.dup_expanded %>% as_tibble() %>%
  dplyr::filter(V1 == V4) %>% # remove chimeric alignments
  rowwise() %>%
  mutate(start = min(V2, V3, V5, V6),
         end = max(V2, V3, V5, V6),
         width = end - start) %>%
  mutate(paired_end_reads = 0,
         mateid = 0,
         mean_mapq = 0) %>%
  dplyr::rename("seqnames"="V1") %>%
  dplyr::select(seqnames, mateid, paired_end_reads, start, end, width, mean_mapq) # %>% head()
# bedpe9_for_medremix %>% head()
# saveRDS(object = bedpe9_for_medremix, file = paste0(output_fpath, ".Rds"))

message('Writing out bedpe4medremix file...')
## { takes a while, ~1hr for 1GB bedpe.gz file }
bedpe9_for_medremix %>%
  write.table(file = output_fpath, sep = "\t", quote = F, row.names = F, col.names = F)

