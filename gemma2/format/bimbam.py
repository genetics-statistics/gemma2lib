# GEMMA2 BIMBAM format support

import json
import gzip
import logging
import numpy as np
from os.path import dirname, basename, splitext
import sys

from types import SimpleNamespace
from gemma2.utility.system import memory_usage

from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno

def write_bimbam(controlfn,compression_level):
    """Write BIMBAM files from R/qtl2 and GEMMA control file"""
    path = dirname(controlfn)
    control = load_control(controlfn)
    base = splitext(control.pheno)[0]
    if path:
        base = path + "/" + base

    phenofn = base+"_bimbam.txt"
    logging.info(f"Writing pheno file {phenofn}")
    with open(phenofn,"w") as f:
        for p in iter_pheno(control.pheno, sep=control.sep, header=False):
            # skip the header and the item counter, otherwise same
            f.write("\t".join(p[1:]))

    base = splitext(splitext(control.geno)[0])[0]
    if path:
        base = path + "/" + base
    genofn = base+"_bimbam.txt.gz"
    logging.info(f"Writing geno file {genofn}")
    genotype_translate = { "A": "1", "B": "0", "H": "2"}
    with gzip.open(genofn, mode='wb', compresslevel=compression_level) as f:
        # f.write("marker".encode())
        for marker,genotypes in iter_geno(control.geno, sep=control.geno_sep, header=False):
            f.write(marker.encode())
            f.write(" - - ".encode())
            f.write(" ".join([genotype_translate[v] for v in genotypes]).encode())
            f.write("\n".encode())
