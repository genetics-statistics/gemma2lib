import logging
import numpy as np
from os.path import dirname, basename
from subprocess import run,CompletedProcess

from gemma2.format.rqtl2 import load_control, load_geno
from gemma2.utility.options import get_options_ns
from gemma2.utility.system import memory_usage

def compute_kinship(control):
    opts = get_options_ns()
    print(control)
    G,markerlist = load_geno(control)
    print(G)
    g = G
    row = g[0]
    memory_usage()
    markers = control.markers
    print(row[0:markers])
    print(row.shape)
    apply maf_filter
    K = np.dot(g,g.T)
    K = K/markers # FIXME, check MAF filter

    if False:
        # If the matrices are in Fortran order then the computations
        # will be faster when using dgemm.  Otherwise, the function
        # will copy the matrix and that takes time.
        from scipy.linalg.blas import dgemm

        gT = g.T
        A = g
        B = gT

        if not A.flags['F_CONTIGUOUS']:
            AA = A.T
            transA = True
        else:
            AA = A
            transA = False

            if not B.flags['F_CONTIGUOUS']:
                BB = B.T
                transB = True
            else:
                BB = B
                transB = False

                K = dgemm(alpha=1.,a=AA,b=BB,trans_a=transA,trans_b=transB)/markers

    print(K)
    memory_usage()
