# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(viridis)

# Set the working directory
setwd("./")

# ---- Section 1: Read and merge methylation and gene expression data ----

# Read methylation data and intron methylation levels
mgbm <- read.table("./list/0.all.accession.txt", header = TRUE)
inib <- read.table("./list/IBM1.intron.methyl.txt")
ins <- inner_join(mgbm, inib, by = c("GSMnum" = "V1"))

# Add combined name column
ins$namecb <- paste(ins$Ecoid, ins$Name, sep = "_")

# Read global methylation levels and join
global <- read.table("./list/summary.addCMT3ex.txt", header = TRUE)
f2 <- inner_join(ins, global, by = c("GSMnum", "Ecoid", "Name"))

# ---- Section 2: Process gene expression data ----

# Open gene expression data file
con <- file("./list/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv", open = "r")
geneex <- c()

# Read data line by line
while ((line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  a <- unlist(strsplit(line, "\\t"))
  if (a[1] == "AT3G07610") {
    geneex <- a
    break
  }
}
close(con)

# Create dataframe from gene expression
com <- data.frame(Ecoid = as.integer(colnames(geneex)[-1]), IBM1ex = as.integer(geneex[-1]))
ins2 <- inner_join(mgbm, com, by = "Ecoid")
ins2$namecb <- paste(ins2$Ecoid, ins2$Name, sep = "_")

# ---- Section 3: Merge and plot methylation and gene expression ----

# Merge and finalize datasets
final <- inner_join(ins2, ins, by = c("GSMnum", "Ecoid", "Name", "CMT3ex", "hchgnum", "namecb"))
f3 <- inner_join(final, ex, by = "Ecoid")

# Plotting the results
p1 <- ggplot(f3, aes(x = scale(V11), y = hchgnum, label = label)) +
  geom_point(size = 3, alpha = 0.5, color = "blue") +
  geom_text() +
  xlab("mCHG and mCG level in the large intron of IBM1") +
  ylab("number of mCHG-gain genes") +
  theme_classic(base_size = 15) +
  scale_color_gradient(low = "yellow", high = "blue")

# Save plots to PDF
pdf('1.intron_methy_isoformratio.v2.pdf', width = 15, height = 6)
plot_grid(p1, labels = 'B', label_size = 12)
dev.off()

# ---- Section 4: Additional plots for isoform quantification results ----

# Read additional data
l1 <- read.table("./list/KakutaniLab.ibm.isoform.txt", header = TRUE)
l2 <- read.table("./list/Sazelab.ibm.isoform.txt", header = TRUE)

# Process and plot data from additional sources
c1 <- l1 %>%
  group_by(Sample) %>%
  summarise(TPM_mean = mean(geneTPM), TPM_sd = sd(geneTPM), long_mean = mean(long),
            long_sd = sd(long), short_mean = mean(short), short_sd = sd(short)) %>%
  pivot_longer(cols = -Sample, names_to = c("type", ".value"), names_sep = "_")

# Plotting TPM and isoform data
p3 <- ggplot(c1, aes(x = Sample, y = mean, fill = type, ymin = mean - sd, ymax = mean + sd)) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  geom_errorbar(width = 0.2, alpha = 1, size = 0.7, position = position_dodge(0.5), color = "black") +
  scale_fill_manual(values = c("#E77577", "#0065A2", "#FFD872"),
                    labels = c("IBM1-L", "IBM1-S", "total")) +
  xlab("KakutaniLab") +
  ylab("IBM1 Gene expression (TPM)") +
  theme_classic(base_size = 15)

# Save the TPM plot to PDF
pdf('2.labtest_isoformratio.pdf', width = 12, height = 6)
plot_grid(p3, label_size = 12)
dev.off()
