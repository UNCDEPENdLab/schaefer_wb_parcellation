source("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/parcel_brain_behavior_functions.R")
setwd("~/Data_Analysis/schaefer_wb_parcellation")
schaefer_7 <- read.csv("Schaefer2018_200Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num)

schaefer_17 <- read.csv("Schaefer2018_200Parcels_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num)

both <- inner_join(schaefer_7, schaefer_17, by="roi_num")
setDT(both)
fill_mask_with_stats(mask_nifti="Schaefer_244_final_2009c_2.3mm.nii.gz", 
                     stat_dt = both, mask_col = "roi_num", 
                     stat_cols = c("net_num7", "net_num17", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer_7_17", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

both %>% filter(network7=="DorsAttn")
both %>% filter(grepl("DorsAttn", network17))


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

#both <- inner_join(schaefer_7, schaefer_17, by="roi_num")
setDT(kong_17)
fill_mask_with_stats(mask_nifti="high_res_originals/Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz", 
                     stat_dt = kong_17, mask_col = "roi_num", 
                     stat_cols = c("net_num17", "roi_num"),
                     stack_along = "stat_cols", img_prefix = "schaefer200_17_kong", out_dir = "~/Downloads/dan_labels",
                     afni_dir = "~/abin", overwrite = T)

