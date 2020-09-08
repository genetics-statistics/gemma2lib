# GEMMA2 R/qtl2 format support

import json
import gzip
import logging
import numpy as np
import sys

from types import SimpleNamespace
from gemma2.utility.system import memory_usage

def load_control(fn: str) -> SimpleNamespace:
    """Load GEMMA2/Rqtl2 style control file"""
    logging.info(f"Reading GEMMA2/Rqtl2 control {fn}")

    data = json.loads(open(fn).read())
    if not "na.strings" in data.values():
        data["na.strings"] = ["NA","nan","-"]

    data["na_strings"] = data["na.strings"]
    logging.info(data)
    control = SimpleNamespace(**data)
    return control

def load_geno(control):
    """GEMMA2/Rqtl2 eager loading of GENO file. Currently only the compact
format is supported

    """
    fn = control.geno
    inds = control.individuals
    markers = control.markers
    genotype_translate = control.genotypes
    assert hasattr(control, 'geno_compact'), "Expect geno_compact set in control file"
    assert markers>inds, f"markers ({markers}) should be larger than individuals ({inds})"
    logging.info(f"Reading GEMMA2/Rqtl2 geno {fn}")
    shape = (markers,inds)
    g = np.empty(shape, dtype=np.float32, order="F")
    print(g.shape)
    # sys.exit(1)
    in_header = True
    with gzip.open(fn) as f:
        line = f.readline()
        if line[0] == '#':
            next
        i = None
        while line:
            if in_header:
                in_header = False
                i = 0
                next
            else:
                (marker,l) = line.decode().rstrip().split("\t",3)
                # print(list(l.rstrip()))
                g[i,:] = [genotype_translate[v] for v in list(l)]
                # print(i,g[i])
                i += 1
            line = f.readline()
    print(g)
    row = g[0]
    memory_usage()
    print(row[0:markers])
    print(row.shape)
    K = np.dot(g,g.T)/markers

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

def iter_pheno(fn: str, sep: str = "\t", header: bool = False):
    """Iter of GEMMA2 pheno file. Returns by line"""
    count = 0
    logging.info(f"Reading GEMMA2/Rqtl2 pheno {fn}")
    with open(fn,"r") as f:
        for line in f:
            count += 1
            if header or count > 1:
                yield line.strip().split(sep)

def iter_geno(fn: str, sep: str = "\t", header: bool = False):
    count = 0
    logging.info(f"Reading GEMMA2/Rqtl2 geno {fn}")
    with gzip.open(fn) as f:
        for line in f:
            count += 1
            if header or count>1:
                l = line.decode()
                marker,genotypes = l.strip().split("\t",2)
                yield marker,[char for char in genotypes]
