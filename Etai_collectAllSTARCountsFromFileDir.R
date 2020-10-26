#collectAllSTARCountsFromFileDir.R

collectAllSTARCountsFromFileDir <- function(workdir="/czlab/etai/ExperimentsData/Stamatis/results_latest") {
  suffix=".ReadsPerGene.out.tab"
  
  files <- list.files(path = workdir, pattern = sprintf("*%s", suffix), full.names = T)
  
  count <- 0
  htmat <- do.call(cbind, lapply(files,
                                 function(f) {
                                   count <- count + 1
                                   message(count, ": ", f)
                                   read.table(f, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)[,1,drop=F] }
  ))
  names(htmat) <- do.call(rbind, lapply(files, function(f) strsplit(basename(f), split = ".gz")[[1]][1]))
  
  htqc <- htmat[grep("^N_", rownames(htmat)), ]
  htmat <- htmat[-grep("^N_", rownames(htmat)), ]
  
  htqc <- rbind(htqc, geneCounts=colSums(htmat))
  return(list(counts=htmat, qc=htqc))
}

                                        
expr1 <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)
qc <- expr1$qc
counts <- expr1$counts
names(counts) <- gsub(".R1.fastq", "", names(counts))
names(qc) <- gsub(".R1.fastq", "", names(qc))

#Generating QC data:
qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)
qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))
rownames(qcbygene) <- c("th1", "th5", "th10")

qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
qc3 <- data.table(qc3)
QC <- qc3
QC <- data.table(sample_id = gsub(".ReadsPerGene.out.tab", "", colnames(qc)), QC)
saveRDS(QC, file = "/singlecellcenter/etai/ExperimentsData/MOUSE_MELANOMA/data/Yao/H7GJWDRXX/QC.rds")

