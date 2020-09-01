# GEMMA2 R/qtl2 format support

import json
import gzip
import logging
import numpy as np

from types import SimpleNamespace

def load_control(fn: str) -> SimpleNamespace:
    """Load GEMMA2/Rqtl2 style control file"""
    data = json.loads(open(fn).read())
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
    assert hasattr(control, 'geno_compact'), "Expect geno_compact set in control file"
    logging.info(f"Reading {fn}")
    shape = (inds,markers)
    g = np.empty(shape, dtype=float, order='C')

    in_header = True
    with gzip.open(fn) as f:
        line = f.readline()
        if line[0] == '#':
            next
        while line:
            if in_header:
                in_header = False
                next
            else:
                a = line.decode().split("\t",3)
                print(a)
            line = f.readline()
