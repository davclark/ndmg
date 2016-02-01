#!/usr/bin/env python

# Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

# run_m2g.py
# Created by Greg Kiar and Will Gray Roncal on 2016-01-27.
# Email: gkiar@jhu.edu, wgr@jhu.edu

from argparse import ArgumentParser
from datetime import datetime
import ndmg.utils as mgu
import ndmg.register as mgr
import ndmg.track as mgt
import ndmg.graph as mgg
import numpy as np
import nibabel as nb


def pipeline(dti, bvals, bvecs, mprage, atlas, mask, labels, aligned_dti,
             tensors, fibers, graph):
    """
    Creates a brain graph from MRI data 
    """
    startTime = datetime.now()
    # Align DTI volumes to Atlas
    print "Aligning volumes..."
    gtab = mgr().dti2atlas(dti, bvals, bvecs, mprage, atlas, aligned_dti)

    print "Beginning tractography..."
    # Compute tensors and track fiber streamlines
    tracks = mgt().eudx_basic(aligned_dti, mask, gtab,
                                       seed_num=1000000)

    print "Generating graph..."
    # Create graphs from streamlines
    labels_im = nb.load(labels)
    graph = mgg(len(np.unique(labels_im.get_data()))-1, labels)
    graph.make_graph(tracks)
    graph.summary()

    print "Saving derivatives..."
    # Save derivatives to disk
    np.savez(fibers, tracks)
    graph.save_graph(graph)

    print "Execution took: " + str(datetime.now() - startTime)
    print "Complete!"
    pass


def main():
    parser = ArgumentParser(description="This is an end-to-end connectome \
                            estimation pipeline from sMRI and DTI images")
    parser.add_argument("dti", action="store", help="Nifti DTI image stack")
    parser.add_argument("bval", action="store", help="DTI scanner b-values")
    parser.add_argument("bvec", action="store", help="DTI scanner b-vectors")
    parser.add_argument("mprage", action="store", help="Nifti T1 MRI image")
    parser.add_argument("atlas", action="store", help="Nifti T1 MRI atlas")
    parser.add_argument("mask", action="store", help="Nifti binary mask of \
                        brain space in the atlas")
    parser.add_argument("labels", action="store", help="Nifti labels of \
                        regions of interest in atlas space")
    parser.add_argument("aligned_dti", action="store", help="Nifti image of\
                        aligned dti volume in atlas space")
    parser.add_argument("tensors", action="store", help="Computed tensors")
    parser.add_argument("fibers", action="store", help="Fiber streamlines")
    parser.add_argument("graph", action="store", help="Produced graphml file")
    result = parser.parse_args()

    pipeline(result.dti, result.bval, result.bvec, result.mprage, result.atlas,
            result.mask, result.labels, result.aligned_dti, result.tensors,
            result.fibers, result.graph)


if __name__ == "__main__":
    main()