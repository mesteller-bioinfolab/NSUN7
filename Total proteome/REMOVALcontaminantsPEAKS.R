REMOVALcontaminantsPEAKS <- function( yourDataset, contaminants = "contaminants.fasta") {
  if (missing(yourDataset)) {
    stop("Please, enter a valid protein/peptide dataset in .csv format")
  }
  if (grepl(".csv", yourDataset) == FALSE) {
    stop("Please, enter a valid protein/peptide dataset in .csv format")
  }
  message("loading your dataset...\n\n")
  dataset <- read.csv(yourDataset, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  cn <- colnames(dataset)
  if (!any(grepl("^Accession$", cn))) {
    stop("Please, rename the column containing the accession codes with the name 'Accession'")
  }
  message("loading the contaminants database...\n\n")
  cont <- readLines(contaminants)
  cont <- cont[grep(">", cont)]
  cont <- substr(cont, 2, 7)
  cont <- cont[1:211]
  message("Removing the contaminants from your Dataset...\n\n")
  p <- dataset[-c(unlist(lapply(cont, function(x) grep(x, dataset$Accession)))), ]
  pt <- paste(gsub(".csv", "", yourDataset), "_ContRemoved.csv...", sep = "")
  message(paste("Writing the file ", pt, "...", sep = ""))
  write.csv(p, file = pt, row.names = FALSE)
}



