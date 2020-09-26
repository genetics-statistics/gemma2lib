import logging
import numpy as np
from typing import List
from os.path import dirname, basename
from subprocess import run,CompletedProcess
import sys

import gemma2.utility.data as data
from gemma2.filters import maf_num_filter
from gemma2.format.rqtl2 import load_control, load_geno
from gemma2.utility.options import get_options_ns
from gemma2.utility.system import memory_usage

# def impute_gs(marker: str, gs: List[float]) -> List[float]:
#     pass

def compute_kinship(control):
    # FIXME: these values are hard coded to develop the algorithm
    miss = 0.05
    maf = 0.01

    def filter_gs_ok(marker: str, gs: List[float]) -> bool:
        # 1. [X] Always apply the MAF filter when reading genotypes
        # 2. [X] Apply missiness filter
        return maf_num_filter(marker,gs,miss,maf)

    opts = get_options_ns()
    print(control)
    G,markerlist = load_geno(control,filter_gs_ok)
    # print(G)
    ctrl = data.methodize(control)
    # print(type(G))

    for idx, gs in enumerate(G):
        mean = np.mean(gs[~np.isnan(gs)]) # skip NAN
        # 3. [ ] Always impute missing data (injecting the row mean) FIXME
        # 4. [X] Always subtract the row mean
        gs = gs - mean
        # 5. [ ] Center the data by row (which is the default option ~-gk 1~)
        G[idx,:] = gs

    markers = ctrl.markers
    K = np.dot(G,G.T)
    # 6. Always scale the matrix dividing by # of SNPs
    K = K/markers
    print(G)
    print(K)
    memory_usage()
