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
import numpy.ndarray as nar


class quality_control(object):

    def __init__(self):
        """
        Enables quality control of fMRI and DTI processing pipelines
        using the ggplot plotting system.
        """
        pass

    def dice_coefficient(a, b):
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

    def mse(imageA, imageB):
        """
        the 'Mean Squared Error' between the two images is the
        sum of the squared difference between the two images;
        NOTE: the two images must have the same dimension
        from http://www.pyimagesearch.com/2014/09/15/python-compare-two-images/
        """
        imageA = imageA/sum(sum(sum(imageA)))
        imageB = imageB/sum(sum(sum(imageB)))
        err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
        err /= float(imageA.shape[0] * imageA.shape[1])

        # return the MSE, the lower the error, the more "similar"
        # the two images are
        return err

    def plot_kdes(self, before, after, steps=1000):
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
        return [kdebef.evaluate(x_grid), kdeaft.evaluate(x_grid)]

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
        The function produces [numtimepoints] plots in the output
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
        slice_ref = refdata[:, :, 100]

        v_bef = np.zeros(timesteps)
        v_aft = np.zeros(timesteps)

        for t in range(0, timesteps):
            v_bef[t] = self.mse(mri_bdata[:, :, :, t], refdata[:, :, :])
            v_aft[t] = self.mse(mri_adata[:, :, :, t], refdata[:, :, :])
            # csv.write(str(v_bef) + "," + str(v_after) + "\n")
            slice_before = mri_bdata[:, :, 100, t]  # TODO EB:replace with mean
            slice_after = mri_adata[:, :, 100, t]

            plt.subplot(1, 2, 1)
            plt.imshow(slice_before, cmap='gray', interpolation='nearest')
            plt.hold(True)
            plt.imshow(slice_ref, cmap='winter', interpolation='nearest',
                       alpha=0.3)
            plt.title("Similarity = %.2E" % v_bef[t])
            plt.subplot(1, 2, 2)
            plt.imshow(slice_after, cmap='gray', interpolation='nearest')
            plt.hold(True)
            plt.imshow(slice_ref, cmap='winter', interpolation='nearest',
                       alpha=0.3)
            plt.title("Similarity = %.2E" % v_aft[t])
            fnamets = outdir + "/" + fname + "_" + str(t) + "_ts"
            plt.savefig(fnamets + ".png")

        kdes = self.compute_kdes(v_bef, v_after)
        plt.plot(kdes[0])
        plt.hold(True)
        plt.plot(kdes[1])
        plt.title(title + " Hellinger Distance = " +
                  str(self.hdist(kdes[0], kdes[1])))
        plt.xlabel('MSE')
        plt.ylabel('Density')
        plt.legend(handles=['before', 'after'])
        plt.savefig(fname + ".png")
        pass
