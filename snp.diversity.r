setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v5.0-refine/3_genetics/")
library("ggplot2")
library(tidyverse)
library(pegas)
library(ape)
library(cowplot)

# this script is used to obtain the SNP density and Tajima's D for all long introns with/without methylation in arabidopsis
# SNP density

df=read.table("./list/snp_outdir/intron.length.snpsites.txt",header=T) %>%
  mutate(snpdensity=nsites/seq_len)

t.test(snpdensity~status,data=df)

p1=ggplot(df,aes(x=status,y=snpdensity,fill=status))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("#e21f26", "#397fb9"),
                    labels=c("methylated introns","unmethylated introns"))+
  xlab(" ")+
  ylab("SNP density")+
  theme_classic(base_size = 15)+
  scale_x_discrete(labels=c("methylated introns","unmethylated introns"))

# Calculate Tajima's D
# Read the FASTA file

# Define the function to calculate Tajima's D and include methylation status
calculate_tajima_d <- function(file_path) {
  # Read the FASTA file
  dna <- read.dna(file_path, format = "fasta")
  
  # Calculate Tajima's D
  tajima_d_result <- tajima.test(dna)
  
  # Determine methylation status based on filename
  if (grepl("\\bmethylated\\b", basename(file_path))) {
    status <- "methylated"
  } else if (grepl("\\bunmethylated\\b", basename(file_path))) {
    status <- "unmethylated"
  } else {
    status <- "unknown"  # Optional: handle cases where neither keyword is found
  }
  
  # Create a data frame with the results
  result_df <- data.frame(
    file = basename(file_path),
    tajima_d = tajima_d_result$D,
    status = status
  )
  
  return(result_df)
}

# Directory containing the SNP files
snp_dir <- "./list/snp_outdir/"

# List all SNP files in the directory
snp_files <- list.files(snp_dir, pattern = "\\.snp\\.list$", full.names = TRUE)

# Apply the function to each SNP file and combine the results into a single data frame
tajima_d_results <- do.call(rbind, lapply(snp_files, calculate_tajima_d))
t.test(tajima_d~status,data=tajima_d_results)

# Print the combined results
print(tajima_d_results)
View(tajima_d_results)
# Save the results to a CSV file (optional)
write.csv(tajima_d_results, file = "tajima_d_results.csv", row.names = FALSE)

p2=ggplot(tajima_d_results,aes(x=status,y=tajima_d,fill=status))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("#e21f26", "#397fb9"),
                    labels=c("methylated introns","unmethylated introns"))+
  xlab(" ")+
  ylab("Tajima'D")+
  theme_classic(base_size = 15)+
  scale_x_discrete(labels=c("methylated introns","unmethylated introns"))

pdf('all.introns.diversity.pdf',
    width=9,
    height=6)
plot_grid(p1, p2,nrow=1)
dev.off()