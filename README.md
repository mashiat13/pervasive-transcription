# pervasive-transcription
```
library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
# BiocManager::install("ggpubr")
# devtools::install_github("kassambara/ggpubr")
#install.packages("tibble")
library(tibble)
library("ggpubr")
# for each intergenic region block, read that into R. See the starting point and the ending point. Use the starting point to divide by 1million. 


readfile <- function(filepath){
  file <- read.table(filepath)
  colnames(file) <- c("chrom","chromStart", "chromEnd", "genomic_class", "coverage")
  file$length <- file$chromEnd - file$chromStart
  file$log2length <- log(file$length, base = 2)
  file$log2coverage <- log(file$coverage, base = 2)
  file$total_coverage <- file$coverage*file$length
  file$log2_total_coverage <- log(file$total_coverage, base = 2)
  return(file)
}
#example path
setwd('/scratch/junzli_root/junzli/mrabbani/adrienne_bulk/RNAseq_Run_14679R_Nov2017/')
path <- 'map_14679X1/split_by_chrom/coverage_by_chrom/intergenic_cov_chr1_X1.bed'

#intergenic
for (i in c(1:7)){
  path <- paste0('map_14679X',i,'/split_by_chrom/coverage_by_chrom/intergenic_cov_chr1_X',i,'.bed')
  file_cov <- readfile(path)
  print(head(file_cov))
  assign(paste0('intergenic_x',i,'_coverage'), file_cov)
}

#exonic
for (i in c(1:7)){
  path <- paste0('map_14679X',i,'/split_by_chrom/coverage_by_chrom/exonic_cov_chr1_X',i,'.bed')
  file_cov <- readfile(path)
  print(head(file_cov))
  assign(paste0('exonic_x',i,'_coverage'), file_cov)
}

#intronic
for (i in c(1:7)){
  path <- paste0('map_14679X',i,'/split_by_chrom/coverage_by_chrom/intronic_cov_chr1_X',i,'.bed')
  file_cov <- readfile(path)
  print(head(file_cov))
  assign(paste0('intronic_x',i,'_coverage'), file_cov)
}


#look at 3 genomeclass from individual data point
head(intergenic_x1_coverage)


combined_intergenic <- cbind(undiff_sptg=intergenic_x1_coverage$log2coverage, diff_sptg1=intergenic_x2_coverage$log2coverage, diff_sptg2=intergenic_x3_coverage$log2coverage, sptc1 = intergenic_x4_coverage$log2coverage, sptc2=intergenic_x5_coverage$log2coverage,spermatid1= intergenic_x6_coverage$log2coverage, spermatid2=intergenic_x7_coverage$log2coverage)

head(combined_intergenic)

combined_exonic <- cbind(undiff_sptg=exonic_x1_coverage$log2coverage, diff_sptg1=exonic_x2_coverage$log2coverage, diff_sptg2=exonic_x3_coverage$log2coverage, sptc1 = exonic_x4_coverage$log2coverage, sptc2=exonic_x5_coverage$log2coverage,spermatid1= exonic_x6_coverage$log2coverage, spermatid2=exonic_x7_coverage$log2coverage)

combined_intronic <- cbind(undiff_sptg=intronic_x1_coverage$log2coverage, diff_sptg1=intronic_x2_coverage$log2coverage, diff_sptg2=intronic_x3_coverage$log2coverage, sptc1 = intronic_x4_coverage$log2coverage, sptc2=intronic_x5_coverage$log2coverage,spermatid1= intronic_x6_coverage$log2coverage, spermatid2=intronic_x7_coverage$log2coverage)


intergenic_cor7 = cor(combined_intergenic, method = "spearman")
dim(intergenic_cor7)
image(intergenic_cor7)

heatmap(intergenic_cor7, main = "intergenic log2coverage correlation chr1", margins = c(8,5), cexRow = 1, cexCol = 1, keep.dendro = FALSE)


exonic_cor7 = cor(combined_exonic, method = "spearman")
dim(exonic_cor7)
image(exonic_cor7)
heatmap(exonic_cor7, main = "Exon log2coverage correlation chr1", margins = c(8,5), cexRow = 1, cexCol = 1, keep.dendro = FALSE)

intronic_cor7 = cor(combined_intronic, method = "spearman")
dim(intronic_cor7)
image(intronic_cor7)
heatmap(intronic_cor7, main = "Intron log2coverage correlation chr1", margins = c(8,5), cexRow = 1, cexCol = 1, keep.dendro = FALSE)

```


![image](https://user-images.githubusercontent.com/54853508/114716173-caed2780-9d01-11eb-9e48-8d6482bb4b1c.png)
![image](https://user-images.githubusercontent.com/54853508/114716219-d6d8e980-9d01-11eb-8050-d42ea352f10f.png)
![image](https://user-images.githubusercontent.com/54853508/114716253-de988e00-9d01-11eb-8480-ed192e0f0ccc.png)
