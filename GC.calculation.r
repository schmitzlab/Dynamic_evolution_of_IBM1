setwd("C:/Users/zywlm/OneDrive - University of Georgia/02_my project/12_ibm1/v5.0-refine/3_genetics/")
library("ggplot2")
library(tidyverse)

# calculate the CG/AT ratio for target sequences

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")

library(Biostrings)

seq_un <- readDNAStringSet("./list/long.unmethylated.intron.fa")
seq_methy <- readDNAStringSet("./list/long.methylated.intron.fa")
#seq_names=names(seq34)
# Initialize a vector to store CG/AT ratios
cg_at_ratios_un <- numeric(length(seq_un))
cg_at_ratios_met <- numeric(length(seq_methy))
# Loop through each sequence
for (i in seq_along(seq_un)) {
  # Get base counts
  counts <- alphabetFrequency(seq_un[i], baseOnly=TRUE)
  
  # Calculate CG/AT ratio
  cg_count <- counts[2] + counts[3]
  at_count <- counts[1] + counts[4]
  
  # Check for zero to avoid division by zero
  if (at_count > 0) {
    cg_at_ratios_un[i] <- cg_count / (at_count+cg_count)
  } else {
    cg_at_ratios_un[i] <- NA  # Assign NA if AT count is zero
  }
}

for (i in seq_along(seq_methy)) {
  # Get base counts
  counts <- alphabetFrequency(seq_methy[i], baseOnly=TRUE)
  
  # Calculate CG/AT ratio
  cg_count <- counts[2] + counts[3]
  at_count <- counts[1] + counts[4]
  
  # Check for zero to avoid division by zero
  if (at_count > 0) {
    cg_at_ratios_met[i] <- cg_count / (at_count+cg_count)
  } else {
    cg_at_ratios_met[i] <- NA  # Assign NA if AT count is zero
  }
}
# Print or return the CG/AT ratios


median(cg_at_ratios_un)

median(cg_at_ratios_met)

t.test(cg_at_ratios_un,cg_at_ratios_met)

seq_arab <- readDNAStringSet("./list/arab.intron.seq.fa")

seq_names=names(seq_arab)
# Initialize a vector to store CG/AT ratios
cg_at_ratios_a <- numeric(length(seq_arab))

# Loop through each sequence
for (i in seq_along(seq_arab)) {
  # Get base counts
  counts <- alphabetFrequency(seq_arab[i], baseOnly=TRUE)
  
  # Calculate CG/AT ratio
  cg_count <- counts[2] + counts[3]
  at_count <- counts[1] + counts[4]
  
  # Check for zero to avoid division by zero
  if (at_count > 0) {
    cg_at_ratios_a[i] <- cg_count / (at_count+cg_count)
  } else {
    cg_at_ratios_a[i] <- NA  # Assign NA if AT count is zero
  }
}

# divided the 34 introns into 2 part based on methylation or not
seq34 <- readDNAStringSet("./list/34.intron.seq.fa")
seq_names=names(seq34)

# Convert the character array to a dataframe
df <- read.table(text = seq_names, header = FALSE, sep = " ", stringsAsFactors = FALSE)

# Assign column names
colnames(df) <- c("intronindex", "geneid", "species", "start", "end")

# Print the dataframe to check the output

info=read.csv("./list/intron.methy.info.csv",header = T)%>%
  separate(col = label, into = c("geneid", "species"), sep = "_(?=[^_]+_[^_]+$)", remove = TRUE, convert = FALSE) %>%
  mutate(status = if_else(status != "no", "methyl", status))

df_com=inner_join(df,info,by=c("geneid","species","start","end")) %>%
  mutate(length=if_else((end-start)>=1000,"long","short"))%>%
  mutate(label=paste(geneid, species,start, end, length,status,sep=" "))

new_name=as.character(df_com$label)
names(seq34)=new_name
sequence_names=names(seq34)
me_seq <- seq34[
  grepl("long", sequence_names, ignore.case = TRUE) & grepl("methyl", sequence_names, ignore.case = TRUE)
]

un_seq=seq34[
  grepl("long", sequence_names, ignore.case = TRUE) & grepl("no", sequence_names, ignore.case = TRUE)
]

cg_at_ratios_un <- numeric(length(un_seq))
cg_at_ratios_met <- numeric(length(me_seq))
# Loop through each sequence
for (i in seq_along(un_seq)) {
  # Get base counts
  counts <- alphabetFrequency(un_seq[i], baseOnly=TRUE)
  
  # Calculate CG/AT ratio
  cg_count <- counts[2] + counts[3]
  at_count <- counts[1] + counts[4]
  
  # Check for zero to avoid division by zero
  if (at_count > 0) {
    cg_at_ratios_un[i] <- cg_count / (at_count+cg_count)
  } else {
    cg_at_ratios_un[i] <- NA  # Assign NA if AT count is zero
  }
}

for (i in seq_along(me_seq)) {
  # Get base counts
  counts <- alphabetFrequency(me_seq[i], baseOnly=TRUE)
  
  # Calculate CG/AT ratio
  cg_count <- counts[2] + counts[3]
  at_count <- counts[1] + counts[4]
  
  # Check for zero to avoid division by zero
  if (at_count > 0) {
    cg_at_ratios_met[i] <- cg_count / (at_count+cg_count)
  } else {
    cg_at_ratios_met[i] <- NA  # Assign NA if AT count is zero
  }
}

mean(cg_at_ratios_met)
mean(cg_at_ratios_un)

t.test(cg_at_ratios_met,cg_at_ratios_un)
