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

# preproc.py
# Created by Eric Bridgeford on 2016-06-20-16.
# Email: ebridge2@jhu.edu

from subprocess import Popen, PIPE
import numpy as np
import nibabel as nb
import sys


class preproc(object):

    def __init__(self):
        """
        Enables preprocessing of single images for single images. Has options
        to perform motion correction (now).
        """
        pass

    def motion_correct(self, mri, corrected_mri, idx):
        """
        Performs motion correction of a stack of 3D images.

        **Positional Arguments:**
            - mri: the 4d (fMRI) image volume as a nifti file.
            - corrected_mri: the corrected and aligned fMRI image volume.
        """
        cmd = "mcflirt -in " + mri + " -out " + corrected_mri +\
            " -plots -refvol " + str(idx)
        print "Executing: " + cmd
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()

    def preprocess(self, mri, preproc_mri, motion_mri, outdir):
        """
        A function to preprocess a stack of 3D images.

        **Positional Arguments:**
            - mri: the 4d (fMRI) image volume as a nifti file.
            - preproc_mri: the 4d (fMRI) preprocessed image volume
                as a nifti image.
            - outdir: the location to place outputs
        """

        mri_name = op.splitext(op.splitext(op.basename(mri))[0])[0]

        s0 = outdir + "/tmp/" + mri_name + "_0slice.nii.gz"
        qc_mc = outdir + "/qc/mc/"

        # TODO EB: decide whether it is advantageous to align to mean image
        self.motion_correct(mri, motion_mri, 0)

        sys.path.insert(0, '..')  # TODO EB: remove this before releasing

        import utils.utils as mgu
        mgu().get_slice(motion_mri, 0, s0)
        from qc.quality_control import quality_control as mgqc
        mgqc.check_alignments(mri, motion_mri, s0, qc_mc, mri_name,
                              title="Motion Correction")

        cmd = "cp " + motion_mri + " " + preproc_mri
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()
        pass
