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

if (length(args) < 2L) stop("Expects 2 images")
if (length(args) == 3L) {
  out_file <- args[3L]
} else {
  out_file <- "overlap.csv"
}

if (length(args) == 4L) {
  out_img2_trimmed <- args[4L]
} else {
  out_img2_trimmed <- "img2_no_overlap.nii.gz"
}


#img1 <- "long_axis_l_split_1mm.nii.gz"
#img2 <- "Schaefer_244_final_1.0mm.nii.gz"

add_nii_ext <- function(str) {
  if (!grepl("\\.nii(\\.gz)?$", str, perl=T)) {
    str <- paste0(str, ".nii")
    if (!file.exists(str)) {
      str <- paste0(str, ".gz")
    }
  }
  return(str)
}

img1 <- add_nii_ext(args[1])
img2 <- add_nii_ext(args[2])


library(RNifti)
library(glue)

stopifnot(file.exists(img1))
stopifnot(file.exists(img2))

x <- readNifti(img1)
y <- readNifti(img2)

stopifnot(identical(dim(x), dim(y)))

xbin <- (1L * (abs(x) > 1e-5)) # convert to 1/0 image
ybin <- (1L * (abs(y) > 1e-5)) # convert to 1/0 image

overlap <- xbin * ybin
nvox_overlap <- sum(overlap)

cat(glue("First image: {img1}\n", .trim = F))
cat(glue("Second image: {img2}\n", .trim = F))
cat(glue("Non-zero voxels in {basename(img1)}: {sum(xbin)}\n", .trim = F))
cat(glue("Non-zero voxels in {basename(img2)}: {sum(ybin)}\n", .trim = F))
cat(glue("Number of voxels that overlap: {nvox_overlap}\n", .trim = F))

if (nvox_overlap > 0L) {
  opos <- which(overlap == 1L, arr.ind=TRUE)
  df_overlap <- data.frame(opos)
  names(df_overlap) <- c("i", "j", "k")
  df_overlap$img1 <- img1
  df_overlap$img2 <- img2
  df_overlap$img1_val <- x[opos]
  df_overlap$img2_val <- y[opos]
  df_overlap <- df_overlap[,c("img1", "img2", "i", "j", "k", "img1_val", "img2_val"), drop=F]
  # print(df_overlap)
  cat(glue("Values in {basename(img1)} that overlap {basename(img2)}:\n", .trim = F))
  print(table(df_overlap$img1_val))
  cat(glue("Values in {basename(img2)} that overlap {basename(img1)}:\n", .trim = F))
  print(table(df_overlap$img2_val))
  write.csv(df_overlap, file = out_file, row.names = F)
  y_no <- y
  y_no[opos] <- 0 # zero out voxels in 2nd image that overlap first
} else {
  y_no <- y
}

writeNifti(y_no, file = out_img2_trimmed)
quit(save = "no", 0, FALSE) # always returns 0 whether there are overlaps or not