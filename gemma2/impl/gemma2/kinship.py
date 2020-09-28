import logging
import math
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
        values = gs[~np.isnan(gs)]
        mean = np.mean(values) # skip NAN

        if idx == 1:
            print("orig",gs)
            print("mean",mean)
        # print(mean,variance)
        # 3. [X] Always impute missing data (injecting the row mean) FIXME
        def f(value: float):
            if np.isnan(value):
                return mean
            return value
        gs = [f(g) for g in gs]
        # 4. [X] Always subtract the row mean serves to "center" the data
        gs -= mean
        if idx == 1:
            print("mean",gs)
        # 5. [X] Center the data by row (which is the default option
        # ~-gk 1~) std is sqrt(var) to normalize each feature value to
        # a z-score.
        # gs /= np.std(gs) # std is always about 1.0 with BXD. assert std != 0.0:
        # genovar = np.sum(gs**2)/len(gs)-(mean**2)
        # genovar = np.var(gs)
        # 0.968858: geno var
        # if genovar != 0:
        # print(genovar)
        # gs /= math.sqrt(genovar)
        genovar = np.var(gs)
        # gs /= math.sqrt(genovar)
        G[idx,:] = gs
        if idx == 1:
            print("genovar",genovar,np.std(gs))
            print("z-scored",gs)

    G = G[:-1, :]
    print("G",G)
    markers = ctrl.markers-1
    K = np.dot(G.T,G)
    print("raw K",K)
    # 6. Always scale the matrix dividing by # of SNPs
    print("scale",markers)
    K /= float(markers-1)
    # print(G)
    print("G dim",G.shape)
    print("K dim",K.shape)
    print(K)
    memory_usage()
