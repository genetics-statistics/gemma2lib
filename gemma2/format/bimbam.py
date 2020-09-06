# GEMMA2 BIMBAM format support

import json
import gzip
import logging
import numpy as np
from os.path import dirname, basename, splitext
import sys

from types import SimpleNamespace
from gemma2.utility.system import memory_usage

from gemma2.format.rqtl2 import load_control, iter_pheno

def write_bimbam(controlfn):
    """Write BIMBAM files from R/qtl2 and GEMMA control file"""
    path = dirname(controlfn)
    control = load_control(controlfn)
    phenofn = splitext(control.pheno)[0]+"_bimbam.txt"
    if path:
        phenofn = path + "/" + phenofn
    logging.info(f"Writing pheno file {phenofn}")
    with open(phenofn,"w") as f:
        for p in iter_pheno(control.pheno, sep=control.sep, header=False):
            f.write("\t".join(p[1:]))
