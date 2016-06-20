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

    def get_brain(self, brain_file):
        """
        Opens a brain data series for a mask, mri image, or atlas.
        Returns a numpy.ndarray representation of a brain.

        **Positional Arguements**
            -brain_file: an object to open the data for a brain.
                         Can be a string (path to a brain file),
                         nibabel.nifti1.nifti1image, or a numpy.ndarray
        """
        if type(brain_file) is np.ndarray:  # if brain passed as matrix
            braindata = brain_file
        else:
            if type(brain_file) is str:  # object is a path
                brain = nb.load(brain_file)
            elif type(brain_file) is nb.nifti1.Nifti1Image:
                brain = brain_file
            else:
                raise TypeError("Mask file is of type " + type(brain_file) +
                                "; accepted types are numpy.ndarray, " +
                                "string, and nibabel.nifti1.Nifti1Image.")
            braindata = brain.get_data()
        return braindata

    def voxel_timeseries(self, fmri_file, mask_file, voxel_file=""):
        """
        Function to extract timeseries for the voxels in the provided
        mask.
        Returns the voxel timeseries as a numpy.ndarray.

        **Positional Arguments**
            - fmri_file: the path to the fmri 4d volume to extract timeseries.
                         can be string, nifti1image, or ndarray
            - mask_file: the path to the binary mask the user wants to extract
                    voxel timeseries over. Can be string, nifti1image, or
                    ndarray
            - voxel_file: the path to the voxel timeseries to be created.
                          Must be string.
        """
        print "Extracting Voxel Timeseries..."

        # load the mask data
        maskdata = self.get_brain(mask_file)
        maskbool = (maskdata > 0)  # extract timeseries for any labelled voxels

        # load the MRI data
        fmridata = self.get_brain(fmri_file)
        voxel_ts = fmridata[maskbool, :]
        if voxel_file:
            np.savez(voxel_file, voxel_ts)
        return voxel_ts

    def roi_timeseries(self, fmri_file, label_file, roits_file=""):
        """
        Function to extract average timeseries for the voxels in each
        roi of the labelled atlas.
        Returns the roi timeseries as a numpy.ndarray.

        **Positional Arguments**
            - fmri_file: the path to the 4d volume to extract timeseries
            - label_file: the path to the labelled atlas containing labels
                    for the voxels in the fmri image
            - roits_file: the path to where the roi timeseries will be saved
        """
        labeldata = self.get_brain(label_file)

        rois = np.sort(np.unique(labeldata[labeldata > 0]))

        fmridata = self.get_brain(fmri_file)

        # initialize so resulting ra is [numrois]x[numtimepoints]
        roi_ts = np.zeros((rois.shape[0], fmridata.shape[3]))

        for roi in rois:
            roi_idx = np.where(rois == roi)[0][0]  # the index of the roi

            roibool = (labeldata == roi)  # get a bool where our voxels in roi
            roi_voxelwisets = fmridata[roibool, :]

            roi_ts[roi_idx, :] = np.mean(roi_voxelwisets, axis=0)

        if roits_file:
            np.savez(roits_file, roi_ts)
        return roi_ts
