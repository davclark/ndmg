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

# quality_control.py
# Created by Eric W Bridgeford on 2016-06-08.
# Email: ebridge2@jhu.edu

import matplotlib.pyplot as plt
from numpy import ndarray as nar
import numpy as np
from scipy.stats import gaussian_kde
from subprocess import Popen, PIPE
import nibabel as nb
import sys
import re


class quality_control(object):

    def __init__(self):
        """
        Enables quality control of fMRI and DTI processing pipelines
        using the ggplot plotting system.
        """
        pass

    def dice_coefficient(self, a, b):
        """
        dice coefficient 2nt/na + nb.
        Code taken from https://en.wikibooks.org/wiki/Algorithm_Implementation
        /Strings/Dice%27s_coefficient#Python
        ** Positional Arguments:**
            - a: the first array.
            - b: the second array.
        """
        a = nar.tolist(nar.flatten(a))
        b = nar.tolist(nar.flatten(b))
        if not len(a) or not len(b):
            return 0.0
        if len(a) == 1:
            a = a + u'.'
        if len(b) == 1:
            b = b + u'.'

        a_bigram_list = []
        for i in range(len(a)-1):
            a_bigram_list.append(a[i:i+2])
        b_bigram_list = []
        for i in range(len(b)-1):
            b_bigram_list.append(b[i:i+2])

        a_bigrams = set(a_bigram_list)
        b_bigrams = set(b_bigram_list)
        overlap = len(a_bigrams & b_bigrams)
        dice_coeff = overlap * 2.0/(len(a_bigrams) + len(b_bigrams))
        return dice_coeff

    def mse(self, imageA, imageB):
        """
        the 'Mean Squared Error' between the two images is the
        sum of the squared difference between the two images;
        NOTE: the two images must have the same dimension
        from http://www.pyimagesearch.com/2014/09/15/python-compare-two-images/
        NOTE: we've normalized the signals by the mean intensity at each
        point to make comparisons btwn fMRI and MPRAGE more viable.
        Otherwise, fMRI signal vastly overpowers MPRAGE signal and our MSE
        would have very low accuracy (fMRI signal differences >>> actual
        differences we are looking for).
        """
        imageA = imageA/np.mean(imageA)
        imageB = imageB/np.mean(imageB)
        err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
        err /= float(imageA.shape[0] * imageA.shape[1])

        # return the MSE, the lower the error, the more "similar"
        # the two images are
        return err

    def compute_kdes(self, before, after, steps=1000):
        """
        A function to plot kdes of the arrays for the similarity between
        the before/reference and after/reference.

        **Positional Arguments:**
            - before: a numpy array before an operational step.
            - after: a numpy array after an operational step.
        """
        x_min = min(min(before), min(after))
        x_max = max(max(before), max(after))
        x_grid = np.linspace(x_min, x_max, steps)
        kdebef = gaussian_kde(before)
        kdeaft = gaussian_kde(after)
        return [x_grid, kdebef.evaluate(x_grid), kdeaft.evaluate(x_grid)]

    def hdist(self, array1, array2):
        """
        A function for computing the hellinger distance between arrays
        1 and 2.

        **Positional Arguments:**
            - array 1: the first array.
            - array 2: the second array.
        """
        return np.sqrt(np.sum((np.sqrt(array1) -
                              np.sqrt(array2)) ** 2)) / np.sqrt(2)

    def check_alignments(self, mri_bname, mri_aname, refname, outdir,
                         fname, title="", bin=False):
        """
        A function for checking alignments between two data items.
        The function produces [numnvols] plots in the output
        directory. Each plot is a multiplot comprised of [numslices]
        rows, and 2 columns. The left column shows the overlap between
        the mri scan and the reference before any operation, and the
        right shows the overlap between the mri scan and the reference
        after an operation has taken place. The title is the DICE
        overlap between the two images.
        An additional plot is produced to create a density estimate
        of the DICE scores for the entire subject before and after
        an operation has taken place. The title is the hellinger distance.

        **Positional Arguments**
            - mri_bname: the 4d mri file before an operation has taken place
            - mri_aname: the 4d mri file after an operation has taken place
            - refname: the 3d file used as a reference for the operation
            - outdir: the output directory for the plots
            - title: a string (such as the operation name) for use in the plot
                     titles (defaults to "")
            - fname: a string to use in the file handle for easy finding
            - bin=binary: a bool indicating whether to binarize the scans
                          to analyze alignment for (defaults to False)
        """
        print "Performing Quality Control for " + title + "..."
        cmd = "mkdir -p " + outdir + "/" + fname + " " +\
              outdir + "/" + fname + "/timecourse"
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        p.communicate()

        # load the nifti images
        mri_before = nb.load(mri_bname)
        mri_after = nb.load(mri_aname)
        reference = nb.load(refname)

        # load the data contained for each image
        mri_bdata = mri_before.get_data()
        mri_adata = mri_after.get_data()
        refdata = reference.get_data()

        timesteps = mri_bdata.shape[3]
        slices = mri_bdata.shape[2]

        # dice_file = outdir + "/" + fname + ".csv"
        # csv = open(dice_file, 'w')
        refslice = round(0.5 * slices)
        slice_ref = refdata[:, :, refslice]

        v_bef = np.zeros(timesteps)
        v_aft = np.zeros(timesteps)

        for t in range(0, timesteps):
            v_bef[t] = self.mse(mri_bdata[:, :, :, t], refdata[:, :, :])
            v_aft[t] = self.mse(mri_adata[:, :, :, t], refdata[:, :, :])
            # csv.write(str(v_bef) + "," + str(v_after) + "\n")
            slice_before = mri_bdata[:, :, refslice, t]
            # TODO EB:replace refslice and time sequence with mean
            slice_after = mri_adata[:, :, refslice, t]

            plt.subplot(1, 2, 1)
            plt.imshow(slice_before, cmap='gray', interpolation='nearest')
            plt.hold(True)
            plt.imshow(slice_ref, cmap='winter', interpolation='nearest',
                       alpha=0.3)
            plt.title("Error = %.2E" % v_bef[t])
            plt.xlabel('Position (mm)')
            plt.ylabel('Position (mm)')
            plt.subplot(1, 2, 2)
            plt.imshow(slice_after, cmap='gray', interpolation='nearest')
            plt.hold(True)
            plt.imshow(slice_ref, cmap='winter', interpolation='nearest',
                       alpha=0.3)
            plt.title("Error = %.2E" % v_aft[t])
            plt.xlabel('Position (mm)')
            plt.ylabel('Position (mm)')
            plt.yaxis.tick_right()
            fnamets = outdir + "/" + fname + "/timecourse/"
            fnamets = fnamets + fname + "_" + str(t) + "_ts"
            plt.savefig(fnamets + ".png")
            plt.clf()

        kdes = self.compute_kdes(v_bef, v_aft)
        fkde = plt.figure()
        axkde = fkde.add_subplot(111)
        axkde.plot(kdes[0], kdes[1])
        axkde.plot(kdes[0], kdes[2])
        axkde.set_title(title + (" Hellinger Distance = %.4E" %
                        self.hdist(kdes[0], kdes[1])))
        axkde.set_xlabel('MSE')
        axkde.set_ylabel('Density')
        axkde.legend(['before, mean = %.2E' % np.mean(v_bef),
                     'after, mean = %.2E' % np.mean(v_aft)])
        fnamekde = outdir + "/" + fname + "_kde"
        fkde.savefig(fnamekde + ".png")
        pass

    def stat_summary(self, mri, mri_raw, mri_mc, mask, title="", fname=""):
        """
        A function for producing a stat summary page, along with
        an image of all the slices of this mri scan.

        **Outputs:**
            - Mean figure: plot for each voxel's mean signal intensity for
              each voxel in the corrected mri's timecourse
            - STdev figure: plot of the standard deviation of each voxel's
              corrected timecourse
            - SNR figure: plot of signal to noise ratio for each voxel in the
              corrected timecourse
            - Slice Wise Intensity figure: averages intensity over slices
              and plots each slice as a line for corrected mri
                - goal: slice mean signal intensity is approximately same
                  throughout all nvols (lines are relatively flat)
            - motion parameters figures (3):
                - 1 figure for rotational, 1 figure for translational motion
                  params, 1 for displacement, calculated with fsl's mcflirt
            - stat summary file:
                - provides some useful statistics for analyzing fMRI data;
                  see resources for each individual statistic for the
                  computation

        **Positional Arguments:**
            - mri: the path to a corrected mri scan to analyze
            - mri_raw: the path to an uncorrected mri scan to analyze
            - title: the title for the file (ie, Registered, Resampled, etc)
            - fname: the name to give the file.
            - tmpname: the name appended to any temporary files created.
        """
        print "Producing Quality Control Summary. \n" +\
            "\tRaw Image: " + mri_raw + "\n" +\
            "\tCorrected Image: " + mri + "\n" +\
            "\tMask: " + mask + "\n"
        mri_im = nb.load(mri)
        mri_dat = mri_im.get_data()
        mri_raw_im = nb.load(mri_raw)

        print "Opened MRI Images."
        sys.path.insert(0, '..')  # TODO EB: remove this before releasing
        from timeseries.timeseries import timeseries as mgtc

        # get voxel timeseries for some computations
        voxel = mgtc().voxel_timeseries(mri_dat, mask)

        # image for mean signal intensity over time
        mri_datmean = np.mean(mri_dat, axis=3)
        fmean = plt.figure()
        mri_datstd = np.std(mri_dat, axis=3)
        fstd = plt.figure()

        # image for slice SNR = mean / stdev
        mri_datsnr = np.divide(mri_datmean, mri_datstd)
        fsnr = plt.figure()

        # image for slice-wise mean intensity
        mri_datmi = np.squeeze(np.apply_over_axes(np.mean,
                                                  mri_dat, (0, 1)))
        fmi = plt.figure()
        axmi = fmi.add_subplot(111)

        depth = mri_dat.shape[2]
        nvols = mri_dat.shape[3]

        nrows = np.ceil(np.sqrt(depth))
        ncols = np.ceil(depth/nrows)

        mri_dat = None  # done with this, so save memory
        # produce figures for each slice in the image
        for d in range(0, depth):
            axmean = fmean.add_subplot(nrows, ncols, d+1)
            axmean.imshow(mri_datmean[:, :, d], cmap='gray',
                          interpolation='nearest')
            axmean.set_xlabel('Position (res)')
            axmean.set_ylabel('Position (res)')
            axmean.set_title('%d slice' % d)

            axstd = fstd.add_subplot(nrows, ncols, d+1)
            axstd.imshow(mri_datstd[:, :, d], cmap='gray',
                         interpolation='nearest')
            axstd.set_xlabel('Position (res)')
            axstd.set_ylabel('Position (res)')
            axstd.set_title('%d slice' % d)

            axsnr = fsnr.add_subplot(nrows, ncols, d+1)
            axsnr.imshow(mri_datsnr[:, :, d], cmap='gray',
                         interpolation='nearest')
            axsnr.set_xlabel('Position (res)')
            axsnr.set_ylabel('Position (res)')
            axsnr.set_title('%d slice' % d)

            axmi.plot(mri_datmi[d, :])

        axmi.set_xlabel('Timepoint')
        axmi.set_ylabel('Mean Intensity')
        axmi.set_title('Mean Slice Intensity')
        axmi.set_xlim((0, nvols))
        fmean.set_size_inches(nrows*8, ncols*8)
        fstd.set_size_inches(nrows*8, ncols*8)
        fsnr.set_size_inches(nrows*8, ncols*8)
        fmean.savefig(fname + "_mean.png")
        fstd.savefig(fname + "_std.png")
        fsnr.savefig(fname + "_snr.png")
        fmi.savefig(fname + "_slice_intens.png")

        par_file = mri_mc + ".par"

        abs_pos = np.zeros((nvols, 6))
        rel_pos = np.zeros((nvols, 6))
        with open(par_file) as f:
            counter = 0
            for line in f:
                abs_pos[counter, :] = [float(i) for i in re.split("\\s+",
                                                                  line)[0:6]]
                if counter > 0:
                    rel_pos[counter, :] = np.subtract(abs_pos[counter, :],
                                                      abs_pos[counter-1, :])
                counter += 1

        trans_abs = np.linalg.norm(abs_pos[:, 3:6], axis=1)
        trans_rel = np.linalg.norm(rel_pos[:, 3:6], axis=1)
        rot_abs = np.linalg.norm(abs_pos[:, 0:3], axis=1)
        rot_rel = np.linalg.norm(rel_pos[:, 0:3], axis=1)

        ftrans = plt.figure()
        axtrans = ftrans.add_subplot(111)
        axtrans.plot(abs_pos[:, 3:6])  # plots the parameters
        axtrans.set_xlabel('Timepoint')
        axtrans.set_ylabel('Translation (mm)')
        axtrans.set_title('Translational Motion Parameters')
        axtrans.legend(['x', 'y', 'z'])
        axtrans.set_xlim((0, nvols))
        ftrans.savefig(fname + "_trans_mc.png")

        frot = plt.figure()
        axrot = frot.add_subplot(111)
        axrot.plot(abs_pos[:, 0:3])
        axrot.set_xlabel('Timepoint')
        axrot.set_ylabel('Rotation (rad)')
        axrot.set_title('Rotational Motion Parameters')
        axrot.legend(['x', 'y', 'z'])
        axrot.set_xlim((0, nvols))
        frot.savefig(fname + "_rot_mc.png")

        fmc = plt.figure()
        axmc = fmc.add_subplot(111)
        axmc.plot(trans_abs)
        axmc.plot(trans_rel)
        axmc.set_xlabel('Timepoint')
        axmc.set_ylabel('Movement (mm)')
        axmc.set_title('Estimated Displacement')
        axmc.legend(['absolute', 'relative'])
        axmc.set_xlim((0, nvols))
        fmc.savefig(fname + "_disp_mc.png")

        fstat = open(fname + "_stat_sum.txt", 'w')
        fstat.write("General Information\n")
        fstat.write("Raw Image Resolution: " +
                    str(mri_raw_im.get_header().get_zooms()[0:3]) + "\n")
        fstat.write("Corrected Image Resolution: " +
                    str(mri_im.get_header().get_zooms()[0:3]) + "\n")
        fstat.write("Number of Volumes: %d" % nvols)

        fstat.write("\n\n")
        fstat.write("Signal  Statistics\n")
        fstat.write("Signal Mean: %.4f\n" % np.mean(voxel))
        fstat.write("Signal Stdev: %.4f\n" % np.std(voxel))
        fstat.write("Number of Voxels: %d\n" % voxel.shape[0])
        fstat.write("Average SNR per voxel: %.4f\n" %
                    np.nanmean(np.divide(np.mean(voxel, axis=1),
                               np.std(voxel, axis=1))))
        fstat.write("\n\n")

        # Motion Statistics
        mean_abs = np.mean(abs_pos, axis=0)  # column wise means per param
        std_abs = np.std(abs_pos, axis=0)
        max_abs = np.max(np.abs(abs_pos), axis=0)
        mean_rel = np.mean(rel_pos, axis=0)
        std_rel = np.std(rel_pos, axis=0)
        max_rel = np.max(np.abs(rel_pos), axis=0)
        fstat.write("Motion Statistics\n")
        fstat.write("Absolute Translational Statistics>>\n")
        fstat.write("Max absolute motion: %.4f\n" % max(trans_abs))
        fstat.write("Mean absolute motion: %.4f\n" % np.mean(trans_abs))
        fstat.write("Number of absolute motions > 1mm: %d\n" %
                    np.sum(trans_abs > 1))
        fstat.write("Number of absolute motions > 5mm: %d\n" %
                    np.sum(trans_abs > 5))
        fstat.write("Mean absolute x motion: %.4f\n" %
                    mean_abs[3])
        fstat.write("Std absolute x position: %.4f\n" %
                    std_abs[3])
        fstat.write("Max absolute x motion: %.4f\n" %
                    max_abs[3])
        fstat.write("Mean absolute y motion: %.4f\n" %
                    mean_abs[4])
        fstat.write("Std absolute y position: %.4f\n" %
                    std_abs[4])
        fstat.write("Max absolute y motion: %.4f\n" %
                    max_abs[4])
        fstat.write("Mean absolute z motion: %.4f\n" %
                    mean_abs[5])
        fstat.write("Std absolute z position: %.4f\n" %
                    std_abs[5])
        fstat.write("Max absolute z motion: %.4f\n" %
                    max_abs[5])

        fstat.write("Relative Translational Statistics>>\n")
        fstat.write("Max relative motion: %.4f\n" % max(trans_rel))
        fstat.write("Mean relative motion: %.4f\n" % np.mean(trans_rel))
        fstat.write("Number of relative motions > 1mm: %d\n" %
                    np.sum(trans_rel > 1))
        fstat.write("Number of relative motions > 5mm: %d\n" %
                    np.sum(trans_rel > 5))
        fstat.write("Mean relative x motion: %.4f\n" %
                    mean_abs[3])
        fstat.write("Std relative x motion: %.4f\n" %
                    std_rel[3])
        fstat.write("Max relative x motion: %.4f\n" %
                    max_rel[3])
        fstat.write("Mean relative y motion: %.4f\n" %
                    mean_abs[4])
        fstat.write("Std relative y motion: %.4f\n" %
                    std_rel[4])
        fstat.write("Max relative y motion: %.4f\n" %
                    max_rel[4])
        fstat.write("Mean relative z motion: %.4f\n" %
                    mean_abs[5])
        fstat.write("Std relative z motion: %.4f\n" %
                    std_rel[5])
        fstat.write("Max relative z motion: %.4f\n" %
                    max_rel[5])

        fstat.write("Absolute Rotational Statistics>>\n")
        fstat.write("Max absolute rotation: %.4f\n" % max(rot_abs))
        fstat.write("Mean absolute rotation: %.4f\n" % np.mean(rot_abs))
        fstat.write("Mean absolute x rotation: %.4f\n" %
                    mean_abs[0])
        fstat.write("Std absolute x rotation: %.4f\n" %
                    std_abs[0])
        fstat.write("Max absolute x rotation: %.4f\n" %
                    max_abs[0])
        fstat.write("Mean absolute y rotation: %.4f\n" %
                    mean_abs[1])
        fstat.write("Std absolute y rotation: %.4f\n" %
                    std_abs[1])
        fstat.write("Max absolute y rotation: %.4f\n" %
                    max_abs[1])
        fstat.write("Mean absolute z rotation: %.4f\n" %
                    mean_abs[2])
        fstat.write("Std absolute z rotation: %.4f\n" %
                    std_abs[2])
        fstat.write("Max absolute z rotation: %.4f\n" %
                    max_abs[2])

        fstat.write("Relative Rotational Statistics>>\n")
        fstat.write("Max relative rotation: %.4f\n" % max(rot_rel))
        fstat.write("Mean relative rotation: %.4f\n" % np.mean(rot_rel))
        fstat.write("Mean relative x rotation: %.4f\n" %
                    mean_rel[0])
        fstat.write("Std relative x rotation: %.4f\n" %
                    std_rel[0])
        fstat.write("Max relative x rotation: %.4f\n" %
                    max_rel[0])
        fstat.write("Mean relative y rotation: %.4f\n" %
                    mean_rel[1])
        fstat.write("Std relative y rotation: %.4f\n" %
                    std_rel[1])
        fstat.write("Max relative y rotation: %.4f\n" %
                    max_rel[1])
        fstat.write("Mean relative z rotation: %.4f\n" %
                    mean_rel[2])
        fstat.write("Std relative z rotation: %.4f\n" %
                    std_rel[2])
        fstat.write("Max relative z rotation: %.4f\n" %
                    max_rel[2])

        fstat.close()
        pass
