#!/usr/bin/env Rscript

#read in command line arguments.
args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed=TRUE, value=TRUE), fixed=TRUE))
argpos <- grep("--args", args, fixed=TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos+1):length(args)]
} else {
  args <- c()
}

schaefer_txt <- args[1L]
library(dplyr)
library(tidyr)

out_file <- sub("(.*)\\..+$", "\\1.csv", schaefer_txt, perl=TRUE)

schaefer_df <- read.table(schaefer_txt)
stopifnot(ncol(schaefer_df) == 6L)
schaefer_df <- schaefer_df[,1:2] # columns 3-6 are just Yeo colors for freesurfer
names(schaefer_df) <- c("roi_num", "label")
schaefer_df <- schaefer_df %>%

  # some labels contain 5 components (anatomical label within network and region). Force this to 4 components
  mutate(label = sub("([^_]+_[^_]+_[^_]+)_([^_]+)_([^_]+)", "\\1_\\2\\3", label, perl=TRUE)) %>%
  separate(col = label, into=c("net7", "hemi", "network", "subregion"), sep = "_") %>%
  mutate(hemi = substr(hemi, 1, 1)) %>%
  select(-net7)

write.csv(schaefer_df, file=out_file, row.names=FALSE)