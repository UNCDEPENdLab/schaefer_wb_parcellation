source("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/parcel_brain_behavior_functions.R")
setwd("~/Data_Analysis/schaefer_wb_parcellation")
#N.B. The roi_num provided by Schaefer in the 7 and 17 order files differs -- the ROI identities/locations have changed
# To circumvent this, we lookup up the CoM of each ROI, then sort the dataset on x,y,z and assign a spatial_roi_num that matches between orders
# The spatial_roi_num then forms the merging/identifying column

library(glue)
library(readr)
library(dplyr)
library(data.table)
nP <- 200

schaefer_7 <- read.csv(glue("labels/Schaefer2018_{nP}Parcels_7Networks_order.csv")) %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num)

# this has the spatial coordinate, spatial_roi_num
schaefer_7_lookup <- read.csv(glue("labels/Schaefer_{nP}_7networks_labels.csv"))

schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
  rename(roi_num7=roi_num, subregion7=subregion)

schaefer_17 <- read.csv(glue("labels/Schaefer2018_{nP}Parcels_17Networks_order.csv")) %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num) %>%
  select(-hemi) # mirrored in 7

# this has the spatial coordinate, spatial_roi_num
schaefer_17_lookup <- read.csv(glue("labels/Schaefer_{nP}_17networks_labels.csv")) %>%
  select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup

schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
  rename(roi_num17=roi_num, subregion17=subregion)

both <- inner_join(schaefer_7, schaefer_17, by="spatial_roi_num") %>%
  select(spatial_roi_num, roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17, everything())

write.csv(both, glue("Schaefer_{nP}_7_17_mapping.csv"), row.names=F)
setDT(both)

both %>% filter(network7 == "DorsAttn" | network17 %in% c("DorsAttnA", "DorsAttnB")) %>% 
  pull(roi_num7) %>%
  diff()

xtabs(~network7 + network17, both)


fill_mask_with_stats(mask_nifti=glue("Schaefer_{nP+44}_final_2009c_2.3mm.nii.gz"), 
                     stat_dt = both, mask_col = "roi_num7", 
                     stat_cols = c("net_num7", "net_num17", "roi_num7", "roi_num17"),
                     stack_along = "stat_cols", img_prefix = glue("schaefer_{nP}_7_17"), out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

both %>% filter(network7=="DorsAttn")
both %>% filter(grepl("DorsAttn", network17))


both %>% group_by(network17) %>% filter(row_number()==1) %>% arrange(net_num17) %>% 
  #select(roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17)
  select(network7, network17, net_num17) #net_num7, 

# 200 versus 400 parcel comparison
schaefer_200 <- read.csv("labels/Schaefer2018_200Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network200=network, net_num200=net_num, subregion200=subregion)

schaefer_400 <- read.csv("labels/Schaefer2018_400Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network400=network, net_num400=net_num, subregion400=subregion)

setDT(schaefer_400)
setDT(schaefer_200)
fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_ants.nii.gz", 
                     stat_dt = schaefer_400, mask_col = "roi_num", 
                     stat_cols = c("net_num400", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer400_7", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm_ants.nii.gz", 
                     stat_dt = schaefer_200, mask_col = "roi_num", 
                     stat_cols = c("net_num200", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer200_7", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)



#### KONG 200

kong_17 <- read.csv("labels/Schaefer2018_200Parcels_Kong2022_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num)

kong_17 %>% group_by(network17) %>% filter(row_number()==1) %>% arrange(net_num17)

setDT(kong_17)
fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz", 
                     stat_dt = kong_17, mask_col = "roi_num", 
                     stat_cols = c("net_num17", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer200_17_kong", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)



# just testing the ability to pass expressions to R_batch_job, rather than character vectors
# ee <- expression(
#   x <- 15,
#   for (i in 1:10) { print(i)}
# )
# 
# x <- expression({
#   hello <- 5
#   for (i in 1:10) {
#     print(i)
#   }
# })
# 
# r_code = quote({
#   f <- rbindlist(lapply(sub_data$file, fread))
#   setkeyv(a_coordinates, c('x', 'y', 'z'))
#   setkeyv(f, c('x', 'y', 'z'))
#   f <- merge(f, a_coordinates)
#   f[, atlas_name := 'Schaefer_400_remap']
#   f[, atlas_value := NULL]
#   out_file <- file.path('/proj/mnhallqlab/users/michael/sceptic_decon/Schaefer_400_remap/deconvolved', paste0('sub', sub_data$sub[1L], '_run', sub_data$run[1L], '_Schaefer_400_remap_2.3mm_deconvolved.csv.gz'))
#   fwrite(f, file=out_file)
# })


dan_mapping <- read.csv("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/dan_200_400_voxel_overlap.csv")

# schaefer labels
s_labels <- read.table("/Users/hallquist/Data_Analysis/schaefer_wb_parcellation/labels/Schaefer2018_400Parcels_7Networks_order.txt") %>% select(1:2) %>%
  setNames(c("roi_400", "network_400")) %>%
  mutate(hemi=if_else(grepl("_LH_", network_400), "L", "R"),
         network_400=sub("7Networks_[LR]H_", "", network_400)) %>%
  extract(col="network_400", into=c("network_400", "schaefer_label_400"), regex="([^_]+)_(.*)")

s_labels_200 <- read.table("/Users/hallquist/Data_Analysis/schaefer_wb_parcellation/labels/Schaefer2018_200Parcels_7Networks_order.txt") %>% select(1:2) %>%
  setNames(c("roi_200", "network_200")) %>%
  mutate(hemi=if_else(grepl("_LH_", network_200), "L", "R"),
         network_200=sub("7Networks_[LR]H_", "", network_200)) %>%
  extract(col="network_200", into=c("network_200", "schaefer_label_200"), regex="([^_]+)_(.*)") %>%
  select(-hemi)


# afni labels
a_labels <- read.csv("/Users/hallquist/Data_Analysis/schaefer_wb_parcellation/labels/Schaefer_400_7networks_labels.csv") %>%
  dplyr::rename(roi_400=roi_num, glasser_label = MNI_Glasser_HCP_v1.0) %>% dplyr::select(roi_400, glasser_label) %>%
  mutate(glasser_label = make.unique(sub("Focus point:\\s+", "", glasser_label, perl=TRUE)))

dan_mapping <- dan_mapping %>%
  left_join(s_labels, by="roi_400") %>%
  left_join(s_labels_200, by="roi_200") %>%
  left_join(a_labels)

write.csv(dan_mapping, file="/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400.csv", row.names=F)
