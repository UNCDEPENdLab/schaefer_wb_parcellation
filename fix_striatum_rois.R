library("oro.nifti")
library("glue")

sz <- Sys.getenv("sz")
if (sz == "") sz <- "2.3"


##Left
x <- readNIfTI(glue("l_striatum_tight_7Networks_{sz}mm.nii.gz"), reorient=FALSE)

x[x==2] <- 1 #2 -> 1
highvals <- x[x > 2]
highvals <- highvals - 2 # 4 ->2; 5 -> 3; etc.
x[x > 2] <- highvals

print(table(x))

writeNIfTI(x, glue("l_striatum_tight_7Networks_{sz}mm"), verbose = TRUE)

##Right
x <- readNIfTI(glue("r_striatum_tight_7Networks_{sz}mm.nii.gz"), reorient=FALSE)
x[x==3] <- 0 #only 9 voxels to work with

#same renumbering as above
x[x==2] <- 1 #2 -> 1
highvals <- x[x > 2]
highvals <- highvals - 2 # 4 ->2; 5 -> 3; etc.
x[x > 2] <- highvals

print(table(x))

writeNIfTI(x, glue("r_striatum_tight_7Networks_{sz}mm"), verbose = TRUE)
