# This script calculates the number of Peptide-spectrum matches (PSMs) for each interactor in every sample.
# The 'input' data are the files "msms.txt' and 'proteinGroups.txt',
# which were obtained by running the MaxQuant software. The MaxQuant output has been deposited in to the proteomeXchange (see PXD038422)

msms <- read.table("msms.txt", header = TRUE, sep = "\t")
proteins <- read.table("proteinGroups.txt", header = TRUE, sep = "\t")

proteins <- proteins[proteins$Potential.contaminant != "+", ]
proteins <- proteins[proteins$Reverse != "+", ]
proteins <- proteins[proteins$Only.identified.by.site != "+", ]

msms <- msms[msms$Reverse != "+", ]

raw.files <- as.character(unique(msms$Raw.file))

## Create 'sample sheet'
# this file is composed by two columns, one with the raw file name and the second with the name of the sample...
# NOTE: This file MUST BE CREATED MANUALLY using the order of the files that appear in the 'raw.files' object!!!!!

samplesheet <- data.frame("raw" = raw.files, "sample name" = c("EV_1", "EV_2", "EV_3", 
                                                               "CCDC9B_1", "CCDC9B_2", "CCDC9B_3"))
#################################################################################################################

samples <- list()
for (i in 1:length(raw.files)) {
  samples[[i]] <- msms[msms$Raw.file == samplesheet$raw[i], ]
}
names(samples) <- samplesheet$sample.name

psm.samples <- list()
for (i in 1:length(samples)) {
  s <- c()
  message(paste("\n\n\nSample ", i, sep = ""))
  for (ii in 1:nrow(proteins)) {
    pb <- txtProgressBar(min = 1, max = nrow(proteins), style = 3)
    p <- samples[[i]]
    p1 <- p[grep(proteins$Majority.protein.IDs[ii], p$Proteins), ]
    s <- c(s, nrow(p1))
    setTxtProgressBar(pb, ii)
  }
  psm.samples[[i]] <- s
}

close(pb)

names(psm.samples) <- names(samples)
proteins.columns <- c("Majority.protein.IDs", "Protein.names", "Gene.names", "Sequence.length", colnames(proteins)[grep("^Peptides.", colnames(proteins))])
proteins.columns.g <- unlist(lapply(proteins.columns, function(x) grep(x, colnames(proteins))))
proteins.1 <- proteins[, proteins.columns.g]

library(xlsx)

for (i in 1:length(samples)) {
  proteins.1 <- cbind.data.frame(proteins.1, psm.samples[[i]])
}

colnames(proteins.1) <- c(colnames(proteins[, proteins.columns.g]), names(psm.samples))
write.xlsx(proteins.1, "CCDC9B_PSM.xlsx", row.names = FALSE)

# NOTE: the 'proteins.1' object will be used by the 'SAINTFileCreator.R' script to generate the input file for the SAINT algorithm.
# SAINT is necessary to determine the statistics of each potential interactor.

################# END OF THE SCRIPT #####################








