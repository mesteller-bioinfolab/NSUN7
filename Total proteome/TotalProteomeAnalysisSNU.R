library(rstudioapi)
library(ggfortify)
library(limma)
library(pca3d)
library(mice)

setwd(dirname(getActiveDocumentContext()$path)) 

# Removing contaminants...

source("REMOVALcontaminantsPEAKS.R")

REMOVALcontaminantsPEAKS("proteins.csv")

# loading protein list...

proteins <- read.csv("proteins_ContRemoved_merged.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# NOTE: The file "proteins_ContRemoved_merged.csv" has been created by adding the column "#Spec TMT16plex" from the file "proteins2.csv" to
#       the file created by the function "REMOVALcontaminantPEAKS.R" called "proteins_ContRemoved.csv". All these files, 'protein2.csv',
#       'proteins_ContRemoved.csv" and 'proteins_ContRemoved_merged.csv' have been added to the PRIDE/proteomeXchange repository (PXD038418).

prot <- proteins[, c(grep("Accession", colnames(proteins)), grep("Description", colnames(proteins)), 
                 grep("Intensity.TMT", colnames(proteins)), grep("#Spec", colnames(proteins)))]

prot <- prot[, -c(15:18)] # Removing unused columns (non-used TMT labels)...

samples <- grep("Intensity.TMT", colnames(prot))

colnames(prot) <- c("Accession", "Description", 
                    "SNU.EV.1", "SNU.EV.2", "SNU.EV.3",
                    "SNU.NSUN7.1", "SNU.NSUN7.2", "SNU.NSUN7.3",
                    "CL.COND1.1", "CL.COND1.2", "CL.COND1.3",
                    "CL.COND2.1", "CL.COND2.2", "CL.COND2.3",
                    "#Spec")

# log2...

prot[, samples][prot[, samples] == 0] <- NA
prot[, samples] <- log2(prot[, samples])
prot[, samples][prot[, samples] < 0] <- NA

## App/diaspp

check <- apply(prot[, samples], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))
ev.check <- apply(prot[, samples[1:3]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))
nsun7.check <- apply(prot[, samples[4:6]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))

app.ev <- prot[ev.check <= 1 & nsun7.check == 3, ]
app.nsun7 <- prot[nsun7.check <= 1 & ev.check == 3, ]

# This code below created an excel file with the proteins that appear or disappear in a given condition...

writexl::write_xlsx(list("Appear in EV" = app.ev,
                         "Appear in NSUN7" = app.nsun7), "App_disapp_proteins.xlsx")


###########################

## imputation

# NOTE: We have imputated only the proteins with more than 75% of valid value in each experimental condition...
#       The algorithm was 'predictive mean matching'

# 75%

snu.ev <- apply(prot[, samples[1:3]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))
snu.nsun7 <- apply(prot[, samples[4:6]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))
cl.cond1 <- apply(prot[, samples[7:9]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))
cl.cond2 <- apply(prot[, samples[10:12]], 1, function(x) length(is.na(x)[is.na(x) == TRUE]))

prot.75 <- prot[snu.ev <= 1 & snu.nsun7 <= 1 & cl.cond1 <= 1 & cl.cond2 <= 1, ]

imp1 <- mice(prot.75[, samples], m = 5, maxit = 50, method = "pmm", seed = 500)

prot.75.imp <- prot.75
imp2 <- complete(imp1, 2)
prot.75.imp[, samples] <- imp2

# Normalization against the median...

prot.filt <- prot.75.imp

medians <- apply(prot.filt[, samples], 2, median, na.rm = TRUE)
median.all <- median(as.matrix(prot.filt[, samples]), na.rm = TRUE)

prot.filt.norm <- t(t(prot.filt[, samples]) - medians) + median.all
prot.filt.norm <- cbind.data.frame(prot.filt[, 1:2], prot.filt.norm)
prot.filt.norm$`#Spec` <- prot.filt[, ncol(prot.filt)]

# PCA

pca4 <- prcomp(t(prot.filt.norm[, samples]), scale = TRUE, center = TRUE)

t_prot.filt.norm <- as.data.frame(t(prot.filt.norm[, samples]))
t_prot.filt.norm$groups <- c(rep("SNU-423-EV", 3), rep("SNU-423-NSUN7", 3), rep("CL-COND1", 3), rep("CL-COND2", 3))

png("PCA_all_samples.png")
autoplot(pca4, data = t_prot.filt.norm, label = TRUE, colour = "groups", frame = TRUE)
dev.off()

# limma

groups <- c(rep("SNU_423_EV", 3), rep("SNU_423_NSUN7", 3), rep("CL_COND1", 3), rep("CL_COND2", 3))

design <- model.matrix(~ 0 + groups)
colnames(design) <- gsub("groups", "", colnames(design))
fit <- lmFit(prot.filt.norm[, samples], design = design)
cnt <- makeContrasts("SNU-423-NSUN7 vs SNU-423-EV" = SNU_423_NSUN7 - SNU_423_EV,
                     levels = design)
fit1 <- contrasts.fit(fit, contrasts = cnt)
fit1 <- eBayes(fit1)

tt <- topTable(fit1, number = nrow(prot.filt.norm))

tt1 <- topTable(fit1, coef = "SNU-423-NSUN7 vs SNU-423-EV", number = nrow(prot.filt.norm))
tt1.idx <- unlist(lapply(rownames(tt1), function(x) which(rownames(prot.filt.norm) == x)))
tt1.f <- cbind.data.frame(prot.filt.norm[tt1.idx, ], tt1)

# 

library(DEqMS)

# NOTE: WE use this package to implement the accuracy of the protein abundance estimation by using the number 
#       of peptides/PSMs quantified for each protein

fit1.1 <- fit1

fit1.1$count <- prot.filt.norm$`#Spec`
fit1.2 <- spectraCounteBayes(fit1.1)

DEqMS.results.snu <- outputResult(fit1.2, coef_col = 1)
snu.idx <- unlist(lapply(rownames(DEqMS.results.snu), function(x) which(rownames(prot.filt.norm) == x)))
DEqMS.results.snu.f <- cbind.data.frame(prot.filt.norm[snu.idx, ], DEqMS.results.snu)

write.csv(DEqMS.results.snu.f, "SNU_DEqMS_merged.csv", row.names = FALSE)

# IMPORTANT NOTE: The file "SNU_DEqMS_merged.csv" has been reordered to obtain the final file "NSUN7_IVNS1ABP_SNUcells.xlsx",
#                 already included in the PRIDE/proteomeXchange repository (PXD038418).

#################### END OF THE SCRIPT ##########################


