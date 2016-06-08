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

    def voxel_timeseries(self, fmri, mask):
        """
        Function to extract timeseries for the voxels in the provided
        mask.

        **Positional Arguments**
            - fmri: the path to the fmri 4d volume to extract timeseries for
            - mask: the path to the binary mask the user wants to extract the
                    voxel timeseries over
        """
        pass

    def roi_timeseries(self, fmri, label):
        """
        Function to extract average timeseries for the voxels in each
        roi of the labelled atlas.

        **Positional Arguments**
            - fmri: the path to the 4d volume to extract timeseries for
            - label: the path to the labelled atlas containing roi labels
                    for the voxels in the fmri image
        """
        pass
