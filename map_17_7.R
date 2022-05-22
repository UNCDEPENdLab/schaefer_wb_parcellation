source("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/parcel_brain_behavior_functions.R")
setwd("~/Data_Analysis/schaefer_wb_parcellation")
#N.B. The roi_num provided by Schaefer in the 7 and 17 order files differs -- the ROI identities/locations have changed
# To circumvent this, we lookup up the CoM of each ROI, then sort the dataset on x,y,z and assign a spatial_roi_num that matches between orders
# The spatial_roi_num then forms the merging/identifying column

schaefer_7 <- read.csv("labels/Schaefer2018_200Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num)

# this has the spatial coordinate, spatial_roi_num
schaefer_7_lookup <- read.csv("labels/Schaefer_200_7networks_labels.csv")

schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
  rename(roi_num7=roi_num, subregion7=subregion)

schaefer_17 <- read.csv("labels/Schaefer2018_200Parcels_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num) %>%
  select(-hemi) # mirrored in 7

# this has the spatial coordinate, spatial_roi_num
schaefer_17_lookup <- read.csv("labels/Schaefer_200_17networks_labels.csv") %>%
  select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup

schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
  rename(roi_num17=roi_num, subregion17=subregion)

both <- inner_join(schaefer_7, schaefer_17, by="spatial_roi_num") %>%
  select(spatial_roi_num, roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17, everything())
setDT(both)


fill_mask_with_stats(mask_nifti="Schaefer_244_final_2009c_2.3mm.nii.gz", 
                     stat_dt = both, mask_col = "roi_num7", 
                     stat_cols = c("net_num7", "net_num17", "roi_num7", "roi_num17"),
                     stack_along = "stat_cols", img_prefix = "schaefer_7_17", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

both %>% filter(network7=="DorsAttn")
both %>% filter(grepl("DorsAttn", network17))


both %>% group_by(network17) %>% filter(row_number()==1) %>% arrange(net_num17) %>% 
  #select(roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17)
  select(network7, network17, net_num17) #net_num7, 

# 400 parcel attempt
schaefer_7 <- read.csv("labels/Schaefer2018_400Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num)

schaefer_17 <- read.csv("labels/Schaefer2018_400Parcels_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num)

both <- inner_join(schaefer_7, schaefer_17, by="roi_num")
setDT(both)
fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_ants.nii.gz", 
                     stat_dt = both, mask_col = "roi_num", 
                     stat_cols = c("net_num7", "net_num17", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer400_7_17", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)



#### KONG 200

kong_17 <- read.csv("labels/Schaefer2018_200Parcels_Kong2022_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num)

kong_17 %>% group_by(network17) %>% filter(row_number()==1) %>% arrange(net_num17)

#both <- inner_join(schaefer_7, schaefer_17, by="roi_num")
setDT(kong_17)
fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz", 
                     stat_dt = kong_17, mask_col = "roi_num", 
                     stat_cols = c("net_num17", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer200_17_kong", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

