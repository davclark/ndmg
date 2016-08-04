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

import os.path as op
from ndmg.utils import utils as mgu
import nibabel as nb
import numpy as np
import nilearn.image as nl
import sys
import dipy.align.reslice as dr
from ndmg.qc import qc as mgqc


class register(object):

    def __init__(self):
        """
        Enables registration of single images to one another as well as volumes
        within multi-volume image stacks. Has options to compute transforms,
        apply transforms, as well as a built-in method for aligning low
        resolution mri images to a high resolution atlas.
        """
        pass

    def align(self, inp, ref, xfm=None, out=None, dof=12, searchrad=True,
              interp="trilinear"):
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
        cmd = "flirt -in " + inp + " -ref " + ref + " -interp " + str(interp)
        if xfm is not None:
            cmd += " -omat " + xfm
        if out is not None:
            cmd += " -out " + out
        if type(dof) is int:
            cmd += " -dof " + str(dof)
        if searchrad:
            cmd += " -searchrx -180 180 -searchry -180 180 " +\
                   "-searchrz -180 180"
        mgu().execute_cmd(cmd)
        pass

    def align_nonlinear(self, inp, ref, xfm, warp, mask=None):
        """
        Aligns two images using nonlinear methods and stores the
        transform between them.

        **Positional Arguments:**
            - inp: the input image.
            - ref: the reference image.
            - affxfm: the affine transform to use.
            - warp: the path to store the nonlinear warp.
        """
        cmd = "fnirt --in=" + inp + " --aff=" + xfm + " --cout=" +\
              warp + " --ref=" + ref + " --subsamp=4,2,1,1"
        if mask is not None:
            cmd += " --refmask=" + mask
        mgu().execute_cmd(cmd)
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
        mgu().execute_cmd(cmd)
        pass

    def apply_warp(self, inp, out, ref, warp, xfm=None, mask=None):
        """
        Applies a warp from the functional to reference space
        in a single step, using information about the structural->ref
        mapping as well as the functional to structural mapping.

        **Positional Arguments:**
            inp: the input image to be aligned as a nifti image file.
            out: the output aligned image.
            ref: the image being aligned to.
            warp: the warp from the structural to reference space.
            premat: the affine transformation from functional to
                structural space.
        """
        cmd = "applywarp --ref=" + ref + " --in=" + inp + " --out=" + out +\
              " --warp=" + warp
        if xfm is not None:
            cmd += " --premat=" + xfm
        if mask is not None:
            cmd += " --mask=" + mask
        mgu().execute_cmd(cmd)
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
        cmd = "eddy_correct " + mri + " " + corrected_mri + " " + str(idx)
        mgu().execute_cmd(cmd)
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
        print "Executing resample for " + base + "to resolution of " +\
              template + "..."
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

    def resample_ant(self, base, res, template):
        """
        A function to resample a base image to that of a template image
        using dipy.
        NOTE: Dipy is far superior for antisotropic -> isotropic
            resampling.

        **Positional Arguments:**
            - base: the path to the base image to resample.
            - res: the filename after resampling.
            - template: the template image to align to.
        """
        print "Resampling..."
        baseimg = nb.load(base)
        tempimg = nb.load(template)
        data2, affine2 = dr.reslice(baseimg.get_data(),
                                    baseimg.get_affine(),
                                    baseimg.get_header().get_zooms()[:3],
                                    tempimg.get_header().get_zooms()[:3])
        img2 = nb.Nifti1Image(data2, affine2)
        print "Saving resampled image..."
        nb.save(img2, res)
        pass

    def resample_fsl(self, base, res, template):
        """
        A function to resample a base image in fsl to that of a template.
        **Positional Arguments:**
            - base: the path to the base image to resample.
            - res: the filename after resampling.
            - template: the template image to align to.
        """
        goal_res = int(nb.load(template).get_header().get_zooms()[0])
        cmd = "flirt -in " + base + " -ref " + template + " -out " +\
              res + " -nosearch -applyisoxfm " + str(goal_res)
        mgu().execute_cmd(cmd)
        pass

    def combine_xfms(self, xfm1, xfm2, xfmout):
        """
        A function to combine two transformations, and output the
        resulting transformation.

        **Positional Arguments**
            - xfm1: the path to the first transformation
            - xfm2: the path to the second transformation
            -xfmout: the path to the output transformation
        """
        cmd = "convert_xfm -omat " + xfmout + " -concat " + xfm1 + " " + xfm2
        mgu().execute_cmd(cmd)

    def mri2atlas(self, mri, mprage, atlas, aligned_mri,
                  aligned_mprage, outdir, opt, qcdir="", **kwargs):
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
            atlas_brain = kwargs["atlas_brain"]
            atlas_mask = kwargs["atlas_mask"]
            s0 = outdir + "/tmp/" + mri_name + "_0slice.nii.gz"
            s0_brain = outdir + "/tmp/" + mri_name + "_0slice_brain.nii.gz"
            mprage_brain = outdir + "/tmp/" + mri_name + "_anat_brain.nii.gz"
            xfm_func2mpr = outdir + "/tmp/" + mri_name + "_xfm_func2mpr.mat"
            xfm_mpr2temp = outdir + "/tmp/" + mri_name + "_xfm_mpr2temp.mat"
            warp_mpr2temp = outdir + "/tmp/" + mri_name +\
                "_warp_mpr2temp.nii.gz"

            sys.path.insert(0, '..')
            mgu().get_slice(mri, 0, s0)  # get the 0 slice and save
            # TODO EB: do we want to align the resampled image?
            mgu().extract_brain(mprage, mprage_brain)
            mgu().extract_brain(s0, s0_brain)
            self.align(s0_brain, mprage_brain, xfm_func2mpr)
            self.align(mprage_brain, atlas_brain, xfm_mpr2temp)
            self.align_nonlinear(mprage, atlas, xfm_mpr2temp,
                                 warp_mpr2temp, mask=atlas_mask)
            self.apply_warp(mri, aligned_mri, atlas, warp_mpr2temp,
                            xfm=xfm_func2mpr)
            self.apply_warp(mprage, aligned_mprage, atlas, warp_mpr2temp,
                            mask=atlas_mask)

            # get a brain mask for the s0 slice and extract the brain
            # mgu().extract_brain(s0, s0_brain, opts="-m")
            # mgu().extract_brain(mprage, mprage_bet)
            # self.align(s0_brain, mprage_bet, xfm_func2mpr)
            # self.align(mprage, atlas, xfm_mpr2temp)
            # self.combine_xfms(xfm_mpr2temp, xfm_func2mpr, xfm_comb)

            # self.applyxfm(mri, atlas, xfm_comb, aligned_mri)
            # self.applyxfm(mprage, atlas, xfm_mpr2temp, aligned_mprage)

            # mgu().apply_mask(aligned_mri, mri_brain, s0_mask)
            # mgu().extract_brain(aligned_mprage, mprage_brain)
            if qcdir is not None:                
                mgqc().check_alignments(mri, aligned_mri, atlas, qcdir,
                                        mri_name, title="Registration")
                mgqc().image_align(aligned_mri, atlas_brain, qcdir,
                                   scanid=mri_name, refid=atlas_name)

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
            self.align_slices(mri, mri2, np.where(gtab.b0s_mask)[0])

            # Loads DTI image in as data and extracts B0 volume
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
