#####################################################################################
###############                 Index Table Generator                 ###############
#####################################################################################
#
# Author: Stefano Secchia
#
# Function: this script generates all the possible cell barcodes (combination of Tn5 and PCR indexes) that are present in a sci-ATAC experiment
# 
# Usage: Rscript indextable_generator.R INPATH SAMPLELIST PLATELIST OUTPATH PREFIX
#
######################################################################################

# Arguments
args = commandArgs(TRUE)

inpath <- args[1]
sample.list <- args[2]
plate.list <- args[3]
outpath <- args[4]
prefix <- args[5]

# Import barcode tables
tn5_i7 <- read.table(paste0(inpath, "/tn5_i7.txt"), header=T, stringsAsFactors=F)
tn5_i7 <- tn5_i7[,c("column", "revComp")]
tn5_i5 <- read.table(paste0(inpath, "/tn5_i5.txt"), header=T, stringsAsFactors=F)
pcr_p7 <- read.table(paste0(inpath, "/pcr_P7.txt"), header=T, stringsAsFactors=F)
pcr_p7 <- pcr_p7[c("row", "revComp")]
pcr_p5 <- read.table(paste0(inpath, "/pcr_P5.txt"), header=T, stringsAsFactors=F)
pcr_p5 <- pcr_p5[,c("column", "index")]

# Generate all possible tn5 barcodes (96 combinations)
tn5_all_barc <- data.frame(matrix(ncol = 4, nrow = 0), stringsAsFactors = F)
colnames(tn5_all_barc) <- c("row","index", "column", "revComp")

for (a in 1:nrow(tn5_i5)) {
  temp <- tn5_i5[a,]
  temp <- temp[rep(seq_len(nrow(temp)), each = 12), ]
  temp <- cbind(temp, tn5_i7)
  tn5_all_barc <- rbind(tn5_all_barc, temp)
  rm(temp)
}
rm(a, tn5_i5, tn5_i7)

colnames(tn5_all_barc) <- c("row_i5", "index_i5", "column_i7", "index_i7")
rownames(tn5_all_barc) <- as.character(1:nrow(tn5_all_barc))
tn5_all_barc$tn5_id <- paste(tn5_all_barc$row_i5, tn5_all_barc$column_i7, sep="")

# Generate all possible PCR barcodes (9216 combinations)
pcr_all_barc <- data.frame(matrix(ncol = 4, nrow = 0), stringsAsFactors = F)
colnames(pcr_all_barc) <- c("row","revComp", "column", "index")

for (a in 1:nrow(pcr_p7)) {
  temp <- pcr_p7[a,]
  temp <- temp[rep(seq_len(nrow(temp)), each = 12), ]
  temp <- cbind(temp, pcr_p5)
  pcr_all_barc <- rbind(pcr_all_barc, temp)
  rm(temp)
}
rm(a, pcr_p5, pcr_p7)

colnames(pcr_all_barc) <- c("row_p7", "index_p7", "column_p5", "index_p5")
pcr_all_barc$pcr_id <- paste(pcr_all_barc$row_p7, pcr_all_barc$column_p5, sep="")

# Import sample and plate lists
sample.list <- readRDS(sample.list)
plate.list <- readRDS(plate.list)

# Generate the combined Tn5 and PCR barcodes used for the experiment
# N.B: the final barcode is constructed with the order: tn5_i7, pcr_p7, pcr_p5, tn5_i5
barcodes <- data.frame(barcode = vector(), name = vector())
for (sample in names(sample.list)) {
  
  tn5_exp_barc <- tn5_all_barc[tn5_all_barc$tn5_id %in% sample.list[[sample]],]
  pcr_exp_barc <- pcr_all_barc[pcr_all_barc$pcr_id %in% unlist(plate.list),]
  
  full_barc <- vector()
  for (k in 1:nrow(tn5_exp_barc)) {
    for (l in 1:nrow(pcr_exp_barc)) {
      full_barc <- c(full_barc, paste(tn5_exp_barc[k,"index_i7"], pcr_exp_barc[l,"index_p7"], pcr_exp_barc[l,"index_p5"], tn5_exp_barc[k,"index_i5"], sep=""))
    }
  }
  
  rm(k, l)
  barcodes <- rbind(barcodes, data.frame(barcode = full_barc, name = rep(sample, length(full_barc))))
  rm(full_barc, tn5_exp_barc, pcr_exp_barc)
  
}
rm(sample)

# Export table of barcodes by sample
table(barcodes$name)

write.table(barcodes, paste0(outpath, "/", prefix, "_sample_indextable.txt"), row.names=F, col.names=F, sep="\t", quote=F)
rm(barcodes)

# Export table of barcodes by plate
barcodes <- data.frame(barcode = vector(), name = vector())
for (plate in names(plate.list)) {
  
  tn5_exp_barc <- tn5_all_barc[tn5_all_barc$tn5_id %in% unlist(sample.list),]
  pcr_exp_barc <- pcr_all_barc[pcr_all_barc$pcr_id %in% plate.list[[plate]],]
  
  full_barc <- vector()
  for (k in 1:nrow(tn5_exp_barc)) {
    for (l in 1:nrow(pcr_exp_barc)) {
      full_barc <- c(full_barc, paste(tn5_exp_barc[k,"index_i7"], pcr_exp_barc[l,"index_p7"], pcr_exp_barc[l,"index_p5"], tn5_exp_barc[k,"index_i5"], sep=""))
    }
  }
  
  rm(k, l)
  barcodes <- rbind(barcodes, data.frame(barcode = full_barc, name = rep(plate, length(full_barc))))
  rm(full_barc, tn5_exp_barc, pcr_exp_barc)
  
}
rm(plate)

table(barcodes$name)

write.table(barcodes, paste0(outpath, "/", prefix, "_plate_indextable.txt"), row.names=F, col.names=F, sep="\t", quote=F)
rm(barcodes)

