#AFNI outputs in RAI coordinates, we need LPI
#system("3dCM -all_rois Schaefer2016_400Parcels_7Networks_colors_23_05_16_fsl152.nii.gz > roi_coords_tmp")
#system("3dCM -all_rois Schaefer2016_400Parcels_7Networks_colors_23_05_16_fsl152.nii.gz > roi_coords_tmp")
#with MH update from march 2018
system("3dCM -all_rois -mask /gpfs/group/mnh5174/default/lab_resources/standard/fsl_mni152/MNI152_T1_2.3mm_brain_mask.nii Schaefer_422_final_jul2018.nii.gz > roi_coords_tmp_422_fsl")

coords <- read.table("roi_coords_tmp_422_fsl")
coords <- coords[-1,] #first row is ROI = 0
names(coords) <- c("x.mni", "y.mni", "z.mni")
coords <- tibble::rownames_to_column(coords, "roi")
coords$roi <- as.numeric(coords$roi) -1 #make numeric and renumber
#str(coords)
#coords$roi <- 1:nrow(coords)
#coords$roi <- paste0("roi", c(1:118, 120:422)) #no 119. match names with adjmats
#coords$roi <- c(1:118, 120:422) #no 119. match names with adjmats

#need to flip X and Y to get RAI -> LPI
coords$x.mni <- -1*coords$x.mni
coords$y.mni <- -1*coords$y.mni

# #need to flip X to get RPI -> LPI
# coords$x.mni <- -1*coords$x.mni
# #coords$y.mni <- -1*coords$y.mni

write.table(coords[,c("roi", "x.mni", "y.mni", "z.mni")], file="Schaefer422_fsl_roi_coords_CoM_lpi.txt", row.names=FALSE, col.names=FALSE)

#lookup the coordinates
system("LookupROINames -in_file Schaefer422_fsl_roi_coords_CoM_lpi.txt -out_file Schaefer422_fsl_roi_coords_CoM_lpi_anatlabels.txt -xyz_cols 2 3 4 -roinum_col 1 -afnidir /opt/aci/sw/afni/17.0.02/bin")

unlink("roi_coords_tmp_422_fsl")
