#!/bin/bash


while read -r rest <&3 && read -r anat <&4; do
       python ndmg/ndmg/fmri_pipeline.py $rest $anat /home/eric/design_team/dependencies/atlas/MNI152_T1_2mm.nii.gz /home/eric/design_team/dependencies/atlas/MNI152_T1_2mm_brain.nii.gz /home/eric/design_team/dependencies/mask/MNI152_T1_2mm_brain_mask.nii.gz /mnt/ssd2/data/NKI-ENH/NKI/outputs/fngs_v714_fnirtv1 /home/eric/design_team/dependencies/label/desikan.nii.gz --fmt graphml
done 3<$1 4<$2
