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
