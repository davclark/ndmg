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

# register.py
# Created by Greg Kiar on 2016-01-28.
# Email: gkiar@jhu.edu, ebridge2@jhu.edu

from subprocess import Popen, PIPE
import os.path as op
import ndmg.utils as ndu
import nibabel as nb
import numpy as np
import nilearn.image as nl
import sys  # remove this before releasing; only here so we can get new utils


class register(object):

    def __init__(self):
        """
        Enables registration of single images to one another as well as volumes
        within multi-volume image stacks. Has options to compute transforms,
        apply transforms, as well as a built-in method for aligning low
        resolution mri images to a high resolution atlas.
        """
        pass

    def align(self, inp, ref, xfm):
        """
        Aligns two images and stores the transform between them

        **Positional Arguments:**

                inp:
                    - Input impage to be aligned as a nifti image file
                ref:
                    - Image being aligned to as a nifti image file
                xfm:
                    - Returned transform between two images
        """
        cmd = "flirt -in " + inp + " -ref " + ref + " -omat " + xfm +\
              " -cost mutualinfo -bins 256 -dof 12 -searchrx -180 180" +\
              " -searchry -180 180 -searchrz -180 180"
        print "Executing: " + cmd
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()
        pass

    def applyxfm(self, inp, ref, xfm, aligned):
        """
        Aligns two images with a given transform

        **Positional Arguments:**

                inp:
                    - Input impage to be aligned as a nifti image file
                ref:
                    - Image being aligned to as a nifti image file
                xfm:
                    - Transform between two images
                aligned:
                    - Aligned output image as a nifti image file
        """
        cmd = "flirt -in " + inp + " -ref " + ref + " -out " + aligned +\
              " -init " + xfm + " -interp trilinear -applyxfm"

        print "Executing: " + cmd
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()
        pass

    def align_slices(self, mri, corrected_mri, idx, opt):
        """
        Performs eddy-correction (or self-alignment) of a stack of 3D images

        **Positional Arguments:**
                mri:
                    - 4D (DTI) image volume as a nifti file
                corrected_mri:
                    - Corrected and aligned DTI volume in a nifti file
                idx:
                    - Index of the volume to align to in the stack. for DTI,
                      this corresponds to the B0 volume.
                opt:
                    - 'f': for fMRI
                    - 'd': for DTI
        """
        if (opt == 'f'):
            cmd = "mcflirt -in " + mri + " -out " + corrected_mri +\
                " -plots -refvol " + str(idx)
            print "Executing: " + cmd
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            p.communicate()
        else:
            cmd = "eddy_correct " + mri + " " + corrected_mri + " " + str(idx)
            print "Executing: " + cmd
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            p.communicate()
        pass

    def resample(self, base, ingested, template):
        """
        Resamples the image such that images which have already been aligned
        in real coordinates also overlap in the image/voxel space.

        **Positional Arguments**
                base:
                    - Image to be aligned
                ingested:
                    - Name of image after alignment
                template:
                    - Image that is the target of the alignment
        """
        # Loads images
        template_im = nb.load(template)
        base_im = nb.load(base)
        # Aligns images
        target_im = nl.resample_img(base_im,
                                    target_affine=template_im.get_affine(),
                                    target_shape=template_im.get_data().shape,
                                    interpolation="nearest")
        # Saves new image
        nb.save(target_im, ingested)
        pass

    def combine_xfms(xfm1, xfm2, xfmout):
        """
        A function to combine two transformations, and output the
        resulting transformation.

        **Positional Arguments**
            - xfm1: the path to the first transformation
            - xfm2: the path to the second transformation
            -xfmout: the path to the output transformation
        """
        cmd = "convert_xfm -omat " + xfmout + " -concat " + xfm1 + " " + xfm2
        print "Executing: " + cmd
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()

    def mri2atlas(self, mri, mprage, atlas, aligned_mri, outdir, opt,
                  **kwargs):
        """
        Aligns two images and stores the transform between them

        **Positional Arguments:**

                mri:
                    - Input impage to be aligned as a nifti image file
                mprage:
                    - Intermediate image being aligned to as a nifti image file
                atlas:
                    - Terminal image being aligned to as a nifti image file
                aligned_mri:
                    - Aligned output mri image as a nifti image file
                opt:
                    -'f' for fMRI
                    -'d' for DTI
                outdir: the base output directory to place files
                **kwargs:
                    -'gtab=tab' the gradient table
                    -'bvals=bval' the bvals
                    -'bvecs=bvec' the bvecs
       """
        # Creates names for all intermediate files used
        # GK TODO: come up with smarter way to create these temp file names
        mri_name = op.splitext(op.splitext(op.basename(mri))[0])[0]
        mprage_name = op.splitext(op.splitext(op.basename(mprage))[0])[0]
        atlas_name = op.splitext(op.splitext(op.basename(atlas))[0])[0]

        if (opt == 'f'):
            mri_mc1 = outdir + "/tmp/" + mri_name + "_mc.nii.gz"
            s0_name = outdir + "/tmp/" + mri_name + "_0slice.nii.gz"
            xfm_0tompr = outdir + "/tmp/" + mri_name + "_xfm_0tompr.mat"
            xfm_mprtotemp = outdir + "/tmp/" + mri_name + "_xfm_mprtotem.mat"
            xfm_comb = outdir + "/tmp/" + mri_name + "_xfm_comb.mat"
            mri_tempreg = outdir + "/tmp/" + mri_name + "_reg.nii.gz"

            # align the fMRI volumes to the 0th volume in each stack
            # EB TODO: figure out whether we want to align to the 0th vol
            # or the mean vol in each stack
            self.align_slices(mri, mri_mc1, 0, 'f')

            mri_mc = nb.load(mri_mc1)

            sys.path.insert(0, '..')  # TODO: remove this before releasing

            import utils.utils as mgu
            s0_im = mgu().get_slice(mri_mc.get_data(), 0)  # get the 0th slice

            # Wraps B0 volume in new nifti image
            s0_head = mri_mc.get_header()
            s0_head.set_data_shape(s0_head.get_data_shape()[0:3])
            s0_out = nb.Nifti1Image(s0_im, affine=mri_mc.get_affine(),
                                    header=s0_head)
            s0_out.update_header()
            nb.save(s0_out, s0_name)

            self.align(s0_name, mprage, xfm_0tompr)
            self.align(mprage, atlas, xfmmprtotemp)
            self.combine_xfms(xfm_mprtotemp, xfm_0tompr, xfm_comb)

            self.applyxfm(mri_mc, atlas, xfm_comb, mri_tempreg)
            self.resample(mri_tempreg, aligned_mri, atlas)

        else:
            gtab = kwargs['gtab']
            bvals = kwargs['bvals']
            bvecs = kwargs['bvecs']

            mri2 = outdir + "/tmp/" + mri_name + "_t2.nii.gz"
            temp_aligned = outdir + "/tmp/" + mri_name + "_ta.nii.gz"
            b0 = outdir + "/tmp/" + mri_name + "_b0.nii.gz"
            xfm1 = outdir + "/tmp/" + mri_name + "_" + mprage_name + "_xfm.mat"
            xfm2 = outdir + "/tmp/" + mprage_name + "_" + atlas_name +\
                "_xfm.mat"
            xfm3 = outdir + "/tmp/" + mri_name + "_" + atlas_name + "_xfm.mat"

            # Align DTI volumes to each other
            self.align_slices(mri, mri2, np.where(gtab.b0s_mask)[0], 'd')

            # Loads DTI image in as data and extracts B0 volume
            import ndmg.utils as mgu
            mri_im = nb.load(mri2)
            b0_im = mgu().get_b0(gtab, mri_im.get_data())
            # GK TODO: why doesn't top import work?

            # Wraps B0 volume in new nifti image
            b0_head = mri_im.get_header()
            b0_head.set_data_shape(b0_head.get_data_shape()[0:3])
            b0_out = nb.Nifti1Image(b0_im, affine=mri_im.get_affine(),
                                    header=b0_head)
            b0_out.update_header()
            nb.save(b0_out, b0)

            # Algins B0 volume to MPRAGE, and MPRAGE to Atlas
            self.align(b0, mprage, xfm1)
            self.align(mprage, atlas, xfm2)

            # Combines transforms from previous registrations in proper order
            self.combine_xfms(xfm2, xfm1, xfm3)

            # Applies combined transform to mri image volume
            self.applyxfm(mri2, atlas, xfm3, temp_aligned)
            self.resample(temp_aligned, aligned_mri, atlas)

            # Clean temp files
            cmd = "rm -f " + mri2 + " " + temp_aligned + " " + b0 + " " +\
                  xfm1 + " " + xfm2 + " " + xfm3
            print "Cleaning temporary registration files..."
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            p.communicate()
