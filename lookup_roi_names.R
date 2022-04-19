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

if (length(args) < 1L) stop("Expects 1 argument: name of nifti containing regions")
if (length(args) == 2L) {
  out_file <- args[2L]
} else {
  out_file <- "roi_labels.csv"
}

library(fmri.pipeline)
library(glue)

#atlas <- "/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_DorsAttn_2.3mm.nii.gz"

atlas <- args[1L]

# these will be in RAI orientation (aka 'DICOM') and the first row is ROI = 0 (not interesting)
system(glue("3dCM -all_rois '{atlas}' > atlas_cm.txt"))

roi_cms <- read.table("atlas_cm.txt")[-1, ] # drop the first row (0 ROI value)
write.table(roi_cms, file = "atlas_cm.txt", row.names = F, col.names = F)

wami_obj <- afni_whereami$new(
  coord_file = "atlas_cm.txt",
  coord_orientation = "RAI", coord_space = "MNI", coord_file_columns = 0:2
)

wami_obj$run()
label_df <- wami_obj$get_whereami_df()
unlink("atlas_cm.txt")
write.csv(label_df, file=out_file, row.names=FALSE)