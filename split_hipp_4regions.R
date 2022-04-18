library(RNifti)

l_hipp <- readNifti("high_res_originals/long_axis_l_1.0mm.nii.gz")
l_hipp_split <- l_hipp

# the cut generates a factor, but when merged back as a numeric, we get
# (posterior) 1 = 0 -- 0.25; 2 = .25-.5; 3 = .5-.75; 4 = .75-1 (anterior)
l_hipp_split[l_hipp_split > 0] <- cut(l_hipp_split[l_hipp_split > 0], breaks=c(0,.25,.5,.75,1))
writeNifti(l_hipp_split, file = "long_axis_l_split_1.0mm.nii.gz")

r_hipp <- readNifti("high_res_originals/long_axis_r_1.0mm.nii.gz")
r_hipp_split <- r_hipp
r_hipp_split[r_hipp_split > 0] <- cut(r_hipp_split[r_hipp_split > 0], breaks=c(0,.25,.5,.75,1))
writeNifti(r_hipp_split, file = "long_axis_r_split_1.0mm.nii.gz")
