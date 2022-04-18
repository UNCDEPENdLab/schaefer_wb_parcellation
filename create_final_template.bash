#!/bin/bash
set -ex

#use the Schaefer 400 variants that were warped to FSL MNI space using the CBIG_Projectfsaverage2MNI_Ants function
#Note that the CBIG_Projectfsaverage2MNI function is listed as obsolete when running CBIG code for Schaefer.
#Furthermore, visual inspection of the parcellation using the old and new projections reveals that the _Ants variants are
#beautifully coregistered, whereas the original variant (from the code on Github) is good but not perfect.
[ ! -d "high_res_originals" ] && mkdir "high_res_originals"

nSchaefer=200 # parcels
nAddl=44 # sum of all subcortical parcels
sz=1.0 # mm: internal size for making combined mask
output_res=2.3 # mm output resolution (resampling)
export sz #picked up in R script
cleanup=1 # whether to remove intermediate files
fsl_to_fonov=1

# 2022: new approach: work in 1.0mm for full combined template, resample at last stage

# N.B. The amount of padding needed in the ResampleImageBySpacing may vary between output resolutions and need manual adjustment!
# Watch out especially for the amygdala parcels.

if [ -n "$1" ]; then
    ref_2009c="$1"
fi

Schaefer_srcdir=/proj/mnhallqlab/lab_resources/parcellation/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
cp ${Schaefer_srcdir}/Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_1mm_ants.nii.gz ./Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_1.0mm_ants.nii.gz

# NB. The Schaefer code produces the parcellation at different resolutions using FreeSurfer's mri_vol2vol. But using ResampleImageBySpacing (from ANTS)
# yields an identical file and can be used throughout this script with 1mm original inputs for consistency.
# Trailing arguments: 0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor
# ResampleImageBySpacing 3 Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_1mm_ants.nii.gz Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_${sz}mm_ants.nii.gz ${sz} ${sz} ${sz} 0 1 1 

#########
#CIT amygdala: 4 ROIs.
#ROI NUMBERING:
# L BLA is nSchaefer + 1
# R BLA is nSchaefer + 2
# L CMN is nSchaefer + 3
# R CMN is nSchaefer + 4

# need to resample CIT to 1mm to align with other 1mm masks
CIT_srcdir=/proj/mnhallqlab/lab_resources/parcellation/CIT168_Amygdala_v1.0.3/CIT168toMNI152_700um

cp ${CIT_srcdir}/CIT168_pAmyNuc_700um_MNI.nii.gz .

#create a BLA mask, which are the sub-briks 1 (BLN_La), 2 (BLN_BL), and 3 (BLN_BM)
fslroi CIT168_pAmyNuc_700um_MNI.nii.gz bln_la 0 1
#fslmaths bln_la -thr 0.4 -bin bln_la -odt char #threshold at 40% and binarize

fslroi CIT168_pAmyNuc_700um_MNI.nii.gz bln_bl 1 1
#fslmaths bln_bl -thr 0.4 -bin bln_bl -odt char #threshold at 40% and binarize

fslroi CIT168_pAmyNuc_700um_MNI.nii.gz bln_bm 2 1
#fslmaths bln_bm -thr 0.4 -bin bln_bm -odt char #threshold at 40% and binarize

#add individual masks and threshold if the voxel has at a summed 40% probability of being the basolateral nucleus
#use -fillh to remove any discontinuity in the mask (internal holes)
fslmaths bln_la -add bln_bl -add bln_bm -thr 0.4 -bin -fillh bla_combined -odt char

imrm bln_la bln_bl bln_bm #remove temporary images for each subnucleus

#separate left and right
3dcalc -overwrite -LPI -a bla_combined.nii.gz -expr 'a*isnegative(x)' -prefix l_bla.nii.gz
3dcalc -overwrite -LPI -a bla_combined.nii.gz -expr 'a*ispositive(x)' -prefix r_bla.nii.gz

#use cen + cmn for mask
fslroi CIT168_pAmyNuc_700um_MNI.nii.gz cen_prob 3 1
#fslmaths cen_prob -thr 0.4 -bin -fillh cen_mask -odt char

#separate left and right
#3dcalc -overwrite -LPI -a cen_mask.nii.gz -expr 'a*isnegative(x)' -prefix l_cen.nii.gz
#3dcalc -overwrite -LPI -a cen_mask.nii.gz -expr 'a*ispositive(x)' -prefix r_cen.nii.gz

#broader centromedial amygdala
fslroi CIT168_pAmyNuc_700um_MNI.nii.gz cmn_prob 4 1
#fslmaths bln_la -thr 0.4 -bin bln_la -odt char #threshold at 40% and binarize

#combine cen and cmn for mask
fslmaths cen_prob -add cmn_prob -thr 0.4 -bin -fillh cmn_combined -odt char

3dcalc -overwrite -LPI -a cmn_combined.nii.gz -expr 'a*isnegative(x)' -prefix l_cmn.nii.gz
3dcalc -overwrite -LPI -a cmn_combined.nii.gz -expr 'a*ispositive(x)' -prefix r_cmn.nii.gz

amy_inputs=(l_bla.nii.gz r_bla.nii.gz l_cmn.nii.gz r_cmn.nii.gz)
for a in "${amy_inputs[@]}"; do
    ResampleImageBySpacing 3 ${a} ${a/.nii.gz/}_${sz}mm.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor
    imrm "$a"
done

fslmaths l_bla_${sz}mm.nii.gz -add ${nSchaefer} -thr $((nSchaefer + 1)) l_bla_${sz}mm_$((nSchaefer + 1))
fslmaths r_bla_${sz}mm.nii.gz -add $((nSchaefer + 1)) -thr $((nSchaefer + 2)) r_bla_${sz}mm_$((nSchaefer + 2))
fslmaths l_cmn_${sz}mm.nii.gz -add $((nSchaefer + 2)) -thr $((nSchaefer + 3)) l_cmn_${sz}mm_$((nSchaefer + 3))
fslmaths r_cmn_${sz}mm.nii.gz -add $((nSchaefer + 3)) -thr $((nSchaefer + 4)) r_cmn_${sz}mm_$((nSchaefer + 4))

immv CIT168_pAmyNuc_700um_MNI high_res_originals

######
#thalamus
## ROIs are X05-X08 for L; X09-X12 for R
Thalamus_srcdir=$FSLDIR/data/atlases/Thalamus
cp ${Thalamus_srcdir}/Thalamus-maxprob-thr50-1mm.nii.gz ./Thalamus-maxprob-thr50-1.0mm.nii.gz

# ResampleImageBySpacing 3 Thalamus-maxprob-thr50-1mm.nii.gz Thalamus-maxprob-thr50-${sz}mm.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor

#ROIs 1, 2, 3 are tiny and should be dropped
fslmaths Thalamus-maxprob-thr50-${sz}mm -thr 3.5 -sub 3 -thr 0 Thalamus-maxprob-thr50-${sz}mm_rm123 #second -thr 0 gets rid of the -3 generated by the -sub command

3dcalc -overwrite -LPI -a Thalamus-maxprob-thr50-${sz}mm_rm123.nii.gz -expr 'a*isnegative(x)' -prefix l_thal_${sz}mm.nii.gz
3dcalc -overwrite -LPI -a Thalamus-maxprob-thr50-${sz}mm_rm123.nii.gz -expr 'a*ispositive(x)' -prefix r_thal_${sz}mm.nii.gz

fslmaths l_thal_${sz}mm -add $((nSchaefer + 4)) -thr $((nSchaefer + 5)) l_thal_${sz}_renumber
fslmaths r_thal_${sz}mm -add $((nSchaefer + 8)) -thr $((nSchaefer + 9)) r_thal_${sz}_renumber

immv Thalamus-maxprob-thr50-1.0mm high_res_originals

#########
#Choi 2012 striatum (per Nate: use tight mask)

#this command from freesurfer will transform the Freesurfer-space 256^3 RSP Choi original to the standard FSL MNI space without interpolation
#the --regheader tells the command to apply an identity transformation such that only the grid extents and orientation are modified to match the --targ (template)
mri_vol2vol --targ $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz --regheader \
	    --mov /proj/mnhallqlab/lab_resources/parcellation/Choi_JNeurophysiol12_MNI152/Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz \
	    --o Choi2012_7Networks_MNI152_1.0mm_TightMask.nii.gz

# ResampleImageBySpacing 3 Choi2012_7Networks_MNI152_1.0mm_TightMask.nii.gz Choi2012_7Networks_MNI152_${sz}mm_TightMask.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor

3dcalc -overwrite -LPI -a Choi2012_7Networks_MNI152_${sz}mm_TightMask.nii.gz -expr 'a*isnegative(x)' -prefix l_striatum_tight_7Networks_${sz}mm.nii.gz
3dcalc -overwrite -LPI -a Choi2012_7Networks_MNI152_${sz}mm_TightMask.nii.gz -expr 'a*ispositive(x)' -prefix r_striatum_tight_7Networks_${sz}mm.nii.gz

#note that ROIs 1 and 3 in the original 1mm space are tiny and largely disappear when resampling to 2.3mm. Thus, drop these explicitly.
# > library(oro.nifti)
# oro.nifti 0.9.1
# > x <- readNIfTI("l_striatum_tight_7Networks_2.3mm.nii.gz")
# > table(x)
# x
# 0           2      4      5      6      7
# 585784    114    164     86    241    265

# > x <- readNIfTI("r_striatum_tight_7Networks_2.3mm.nii.gz")
# > table(x)
# x
# 0           2      3      4      5      6      7
# 585834     97      9    189     80    306    139

#this is accomplished using this little fixer script
R CMD BATCH --no-save --no-restore fix_striatum_rois.R

# Choi readout (after dropping 1=Visual and 3=DAN)
# 1 = Somatomotor
# 2 = Salience/VentAttn
# 3 = Limbic
# 4 = FP Control
# 5 = Default

#result: 4 l striatum ROIs are X13-X17; 4 r striatum ROIs are X18-X22
fslmaths l_striatum_tight_7Networks_${sz}mm -add $((nSchaefer + 12)) -thr $((nSchaefer + 13)) l_striatum_tight_7Networks_${sz}mm_renumber
fslmaths r_striatum_tight_7Networks_${sz}mm -add $((nSchaefer + 17)) -thr $((nSchaefer + 18)) r_striatum_tight_7Networks_${sz}mm_renumber

imrm l_striatum_tight_7Networks_${sz}mm r_striatum_tight_7Networks_${sz}mm
immv Choi2012_7Networks_MNI152_1.0mm_TightMask.nii.gz high_res_originals

#verified that inputs match template for ALL final files. Example.
#3dinfo -header_line -same_all_grid Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_2.3mm_ants.nii.gz /gpfs/group/mnh5174/default/lab_resources/standard/fsl_mni152/MNI152_T1_2.3mm.nii

# hippocampus -- break up into posterior tail, posterior body, anterior body, anterior head
# using split along long axis into 4 regions. These are in MNI 2009c space already and
# calculated using the get_hipp_voxels.R script from the Dombrovski et al., 2020 hipp paper.

# create the long axis mask in FSL 1mm space
R CMD BATCH --no-save --no-restore get_hippo_voxels.R

# split each hemisphere into 4 hippocampal parcels
R CMD BATCH --no-save --no-restore split_hipp_4regions.R

# ResampleImageBySpacing 3 long_axis_l_split_1.0mm.nii.gz long_axis_l_split_${sz}mm.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor
# ResampleImageBySpacing 3 long_axis_r_split_1.0mm.nii.gz long_axis_r_split_${sz}mm.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor

# result: L hippo will be X23-X26; R hippo will be X27-X30
fslmaths long_axis_l_split_${sz}mm -add $((nSchaefer + 22)) -thr $((nSchaefer + 23)) long_axis_l_split_${sz}mm_renumber
fslmaths long_axis_r_split_${sz}mm -add $((nSchaefer + 26)) -thr $((nSchaefer + 27)) long_axis_r_split_${sz}mm_renumber

#fslmaths l_striatum_tight_7Networks_${sz}mm -add $((nSchaefer + 12)) -thr $((nSchaefer + 13)) l_striatum_tight_7Networks_${sz}mm_renumber

# add cerebellum

# get buckner atlas to match fonov (actually, do this altogether at the end)
# 3dresample -input atl-Buckner7_space-MNI_dseg.nii -prefix atl-Buckner7_space-MNI_dseg_fonov_1.0mm.nii -master /proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c.nii

buckner7=/proj/mnhallqlab/lab_resources/parcellation/cerebellar_atlases/Buckner_2011/atl-Buckner7_space-MNI_dseg.nii

#this command from freesurfer will transform the Freesurfer-space 256^3 RSP Buckner original to the standard FSL MNI space without interpolation
#the --regheader tells the command to apply an identity transformation such that only the grid extents and orientation are modified to match the --targ (template)
mri_vol2vol --targ $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz --regheader \
	    --mov high_res_originals/atl-Buckner7_space-MNI_dseg.nii.gz  \
	    --o buckner7_fslmni_1.0mm.nii.gz

# 7-network naming is the same as other files
# 1 = visual
# 2 = somatomotor
# 3 = DAN
# 4 = salience/ vent attn
# 5 = limbic
# 6 = FP control
# 7 = default

# ResampleImageBySpacing 3 buckner7_fslmni.nii_1.0mm.gz buckner7_fslmni_${sz}mm.nii.gz ${sz} ${sz} ${sz} 0 1 1 #0 no smoothing, 1 for 1 voxel padding, 1 for nearest neighbor

# divide into l/r
3dcalc -overwrite -LPI -a buckner7_fslmni_${sz}mm.nii.gz -expr 'a*isnegative(x)' -prefix l_buckner7_fslmni_${sz}mm.nii.gz
3dcalc -overwrite -LPI -a buckner7_fslmni_${sz}mm.nii.gz -expr 'a*ispositive(x)' -prefix r_buckner7_fslmni_${sz}mm.nii.gz

# cerebellum will be L X31-X37; R will be X38-X44
fslmaths l_buckner7_fslmni_${sz}mm -add $((nSchaefer + 30))  -thr $((nSchaefer + 31)) l_buckner7_fslmni_${sz}mm
fslmaths r_buckner7_fslmni_${sz}mm -add $((nSchaefer + 37))  -thr $((nSchaefer + 38)) r_buckner7_fslmni_${sz}mm

## COMBINE ALL MASKS

# order matters -- earlier mask values are preserved when there are overlaps
imgs=( l_bla_${sz}mm_$((nSchaefer + 1)) r_bla_${sz}mm_$((nSchaefer + 2)) \
	     l_cmn_${sz}mm_$((nSchaefer + 3)) r_cmn_${sz}mm_$((nSchaefer + 4)) \
	     l_thal_${sz}_renumber r_thal_${sz}_renumber \
	     l_striatum_tight_7Networks_${sz}mm_renumber r_striatum_tight_7Networks_${sz}mm_renumber \
	     long_axis_l_split_${sz}mm_renumber long_axis_r_split_${sz}mm_renumber \
	     l_buckner7_fslmni_${sz}mm r_buckner7_fslmni_${sz}mm \
	     Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_${sz}mm_ants )

fslmaths ${imgs[0]} compile0 # start with base image
for (( i=1; i<${#imgs[@]}; i++ ));
do
    ./check_overlap.R compile$((i-1)).nii.gz ${imgs[$i]} overlap${i}.csv || : #ignore non-zero exit status
    echo "./check_overlap.R compile$((i-1)).nii.gz ${imgs[$i]} overlap${i}.csv"
    ex=$?
    if [ $ex -eq 1 ]; then
      echo "compile failed. Check overlap${i}.csv"
    fi

    #fslmaths compile$((i-1)) -add ${imgs[$i]} compile$i

    # add the non-overlapping component of the 2nd image (give precedence to earlier image)
    fslmaths compile$((i-1)) -add img2_no_overlap compile$i
done

imcp compile$((i-1)) Schaefer_$((nSchaefer + 44))_final_${sz}mm

# can produce unexpected mask values when there are overlaps
# fslmaths Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_${sz}mm_ants \
# 	 -add l_bla_${sz}mm_$((nSchaefer + 1)) -add r_bla_${sz}mm_$((nSchaefer + 2)) \
# 	 -add l_cmn_${sz}mm_$((nSchaefer + 3)) -add r_cmn_${sz}mm_$((nSchaefer + 4)) \
# 	 -add l_thal_${sz}_renumber -add r_thal_${sz}_renumber \
# 	 -add l_striatum_tight_7Networks_${sz}mm_renumber -add r_striatum_tight_7Networks_${sz}mm_renumber \
# 	 -add long_axis_l_split_${sz}mm_renumber -add long_axis_r_split_${sz}mm_renumber \
# 	 -add l_buckner7_fslmni_${sz}mm -add r_buckner7_fslmni_${sz}mm \
# 	 Schaefer_$((nSchaefer + 44))_final_${sz}mm

# resample to final output resolution using voxel size
# in some cases, we have to 3dResample to RPI, then ResampleImageBySpacing the 3dResample to LPI to get the center of the grid to match
ResampleImageBySpacing 3 Schaefer_$((nSchaefer + 44))_final_${sz}mm.nii.gz Schaefer_$((nSchaefer + 44))_final_${output_res}mm.nii.gz ${output_res} ${output_res} ${output_res} 0 1 1 
#3dresample -overwrite -input Schaefer_$((nSchaefer + 44))_final_${sz}mm.nii.gz \
#  -prefix Schaefer_$((nSchaefer + 44))_final_${output_res}mm.nii.gz \
#  -rmode NN -dxyz  ${output_res} ${output_res} ${output_res} 
#  #-master /proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii # (not very general)

# support arbitrary output to 2009c space using a reference image
if [ -n "$ref_2009c" ]; then 
    applywarp -i Schaefer_$((nSchaefer + nAddl))_final_${sz}mm -o Schaefer_$((nSchaefer + nAddl))_final_2009c_fromref --interp=nn \
	      -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef -r "$ref_2009c"
else
    applywarp -i Schaefer_$((nSchaefer + nAddl))_final_${sz}mm -o Schaefer_$((nSchaefer + nAddl))_final_2009c_${sz}mm --interp=nn \
	  -w $MRI_STDDIR/fsl_mni152/fsl_mni152_to_fonov_mni152_warpcoef -r $MRI_STDDIR/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_${sz}mm
fi


#sanity check on final file
#library(oro.nifti)
#x <- readNIfTI("Schaefer_422_final_jul2018.nii.gz", reorient=F)
#sort(unique(as.vector(x)))
#all.equal(sort(unique(as.vector(x))), 0:422)
#TRUE

#Thus, no overlaps leading to ROIs outside the 1-422 range
#CHECK ME IN DETAIL:
#fsleyes $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz $MRI_STDDIR/fsl_mni152/MNI152_T1_2.3mm.nii Schaefer_422_final_jul2018.nii.gz high_res_originals/*

if [ $cleanup -eq 1 ]; then
    imrm l_bla_${sz}mm r_bla_${sz}mm l_cmn_${sz}mm r_cmn_${sz}mm bla_combined cmn_combined cen_prob cmn_prob
    imrm l_bla_${sz}mm_$((nSchaefer + 1)) r_bla_${sz}mm_$((nSchaefer + 2)) l_cmn_${sz}mm_$((nSchaefer + 3)) r_cmn_${sz}mm_$((nSchaefer + 4))
    imrm l_thal_${sz}mm r_thal_${sz}mm Thalamus-maxprob-thr50-${sz}mm_rm123 Thalamus-maxprob-thr50-${sz}mm
    imrm l_thal_${sz}_renumber r_thal_${sz}_renumber
    imrm l_striatum_tight_7Networks_${sz}mm_renumber r_striatum_tight_7Networks_${sz}mm_renumber 
    imrm Choi2012_7Networks_MNI152_${sz}mm_TightMask
    imrm Schaefer2018_${nSchaefer}Parcels_7Networks_order_FSLMNI152_${sz}mm_ants 
    imrm long_axis_l_split_${sz}mm_renumber long_axis_r_split_${sz}mm_renumber
    imrm long_axis_l_split_${sz}mm long_axis_r_split_${sz}mm
    imrm long_axis_l_split_1.0mm long_axis_r_split_1.0mm
    imrm l_buckner7_fslmni_${sz}mm r_buckner7_fslmni_${sz}mm buckner7_fslmni_${sz}mm
    imrm img2_no_overlap
    imrm compile*
fi

cat << EOF > region_labels.txt
# Region labels
1 - ${nSchaefer}: Schaefer cortical parcels
$((nSchaefer + 1)): L BLA
$((nSchaefer + 2)): R BLA
$((nSchaefer + 3)): L CMN
$((nSchaefer + 4)): R BLA
$((nSchaefer + 5)): L Thalamus (need to clarify)
$((nSchaefer + 6)): L Thalamus (need to clarify)
$((nSchaefer + 7)): L Thalamus (need to clarify)
$((nSchaefer + 8)): L Thalamus (need to clarify)
$((nSchaefer + 9)): R Thalamus (need to clarify)
$((nSchaefer + 10)): R Thalamus (need to clarify)
$((nSchaefer + 11)): R Thalamus (need to clarify)
$((nSchaefer + 12)): R Thalamus (need to clarify)
$((nSchaefer + 13)): L Striatum: somatomotor
$((nSchaefer + 14)): L Striatum: salience/vent attn
$((nSchaefer + 15)): L Striatum: limbic
$((nSchaefer + 16)): L Striatum: fp control
$((nSchaefer + 17)): L Striatum: default
$((nSchaefer + 18)): R Striatum: somatomotor
$((nSchaefer + 19)): R Striatum: salience/vent attn
$((nSchaefer + 20)): R Striatum: limbic
$((nSchaefer + 21)): R Striatum: fp control
$((nSchaefer + 22)): R Striatum: default
$((nSchaefer + 23)): L Hippocampus: posterior tail
$((nSchaefer + 24)): L Hippocampus: posterior body
$((nSchaefer + 25)): L Hippocampus: anterior body
$((nSchaefer + 26)): L Hippocampus: anterior head
$((nSchaefer + 27)): R Hippocampus: posterior tail
$((nSchaefer + 28)): R Hippocampus: posterior body
$((nSchaefer + 29)): R Hippocampus: anterior body
$((nSchaefer + 30)): R Hippocampus: anterior head
$((nSchaefer + 31)): L cerebellum: visual
$((nSchaefer + 32)): L cerebellum: somatomotor
$((nSchaefer + 33)): L cerebellum: dan
$((nSchaefer + 34)): L cerebellum: salience/vent attn
$((nSchaefer + 35)): L cerebellum: limbic
$((nSchaefer + 36)): L cerebellum: fp control
$((nSchaefer + 37)): L cerebellum: default
$((nSchaefer + 38)): R cerebellum: visual
$((nSchaefer + 39)): R cerebellum: somatomotor
$((nSchaefer + 40)): R cerebellum: dan
$((nSchaefer + 41)): R cerebellum: salience/vent attn
$((nSchaefer + 42)): R cerebellum: limbic
$((nSchaefer + 43)): R cerebellum: fp control
$((nSchaefer + 44)): R cerebellum: default

The current working directory is: $PWD
You are logged in as: $(whoami)
EOF
