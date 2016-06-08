opyright 2016 NeuroData (http://neurodata.io)
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


import ggplot as gg

class quality_control(object):

    def __init__(self):
        """
        Enables quality control of fMRI and DTI processing pipelines
        using the ggplot plotting system.
        """
        pass

    def check_alignments(self, mri_before, mri_after, reference, outdir,
                        title="", bin=False):
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
            - mri_before: the mri scan before an operation has taken place
            - mri_after: the mri scan after an operation has taken place
            - reference: the scan to use as a reference for the operation
            - outdir: the output directory for the plots
            - title: a string (such as the operation name) for use in the plot
                     titles (defaults to "")
            - bin=binary: a bool indicating whether to binarize the scans
                          to analyze alignment for (defaults to False)
        """
        pass
