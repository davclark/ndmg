#!/usr/bin/env python

# Copyright 2016 NeuroData (http://neurodata.io)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# timeseries.py
# Created by Eric W Bridgeford on 2016-06-07.
# Email: ericwb95@gmail.com

import numpy as np
import nibabel as nb


class timeseries():

    def __init__(self):
        """
        Timeseries extraction class
        """
        pass

    def get_mask(self, mask_ra, value):
        """
        Function to create a brain mask by thresholding the labels

        **Positional Arguements**
            -mask_ra: the mask array file containing the array data
            -value: the value to get the mask for
        """
        pass

    def voxel_timeseries(self, fmri_file, mask_file, voxel_file):
        """
        Function to extract timeseries for the voxels in the provided
        mask.

        **Positional Arguments**
            - fmri_file: the path to the fmri 4d volume to extract timeseries
            - mask_file: the path to the binary mask the user wants to extract
                    voxel timeseries over
            - voxel_file: the path to the voxel timeseries to be created
        """
        print "Extracting Voxel Timeseries..."
        mask = nb.load(mask_file)
        maskdata = mask.get_data()

        maskbool = (maskdata > 0)  # extract timeseries for any labelled voxels

        fmri = nb.load(fmri_file)
        fmridata = fmri.get_data()

        voxel_ts = fmridata[maskbool, :]
        np.savez(voxel_file, voxel_ts)
        pass

    def roi_timeseries(self, fmri_file, label_file, roits_file):
        """
        Function to extract average timeseries for the voxels in each
        roi of the labelled atlas.

        **Positional Arguments**
            - fmri_file: the path to the 4d volume to extract timeseries
            - label_file: the path to the labelled atlas containing labels
                    for the voxels in the fmri image
            - roits_file: the path to where the roi timeseries will be saved
        """
        label = nb.load(label_file)
        labeldata = label.get_data()

        rois = np.sort(np.unique(labeldata[labeldata > 0]))

        fmri = nb.load(fmri_file)
        fmridata = fmri.get_data()

        # initialize so resulting ra is [numrois]x[numtimepoints]
        roi_ts = np.zeros((rois.shape[0], fmridata.shape[3]))

        for roi in rois:
            roi_idx = np.where(rois == roi)[0][0]  # the index of the roi

            roibool = (labeldata == roi)  # get a bool where our voxels in roi
            roi_voxelwisets = fmridata[roibool, :]

            roi_ts[roi_idx, :] = np.mean(roi_voxelwisets, axis=0)

        np.savez(roits_file, roi_ts)
        return roi_ts
