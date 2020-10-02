# GEMMA2 R/qtl2 format support

import json
import gzip
import logging
import numpy as np
import sys
from typing import Callable

from os.path import dirname, basename, isfile
from types import SimpleNamespace
from gemma2.utility.data import methodize
from gemma2.utility.options import get_options_ns
from gemma2.utility.system import memory_usage
import gemma2.utility.safe as safe

def load_control(fn: str) -> dict:
    """Load GEMMA2/Rqtl2 style control file"""
    logging.info(f"Reading GEMMA2/Rqtl2 control {fn}")

    data = json.loads(open(fn).read())
    if not "na.strings" in data.values():
        data["na.strings"] = ["NA","nan","-"]

    data["na_strings"] = data["na.strings"] # does not work as a method
    data["name"] = fn
    logging.info(data)
    return data


def write_new_control(control: dict, transformation: dict):
    """Writes a new control reusing the old version, but adding the last
    command

    """
    opts = get_options_ns()
    cmd = " ".join(opts.args)
    transformation['command'] = cmd
    control['transformations'].append(transformation)
    with safe.control_write_open() as controlf:
        json.dump(control, controlf, indent=4)

def write_control(descr:str, inds:int, markers:int, phenotypes:int, genofn:str,
                  phenofn:str, gmapfn:str, transformation: dict):
    gnfn = basename(genofn)
    phfn = basename(phenofn)
    gmapfn = basename(gmapfn)
    Null = None
    control = {
        "command": descr,
        "crosstype": Null,   # we are not assuming a cross for GEMMA
        "sep": "\t",
        "na.strings": ["-"],
        "comment.char": "#",
        "individuals": inds,
        "markers": markers,
        "phenotypes": phenotypes,
        "geno": gnfn,
        "pheno": phfn,
        "gmap": gmapfn,
        "alleles": ["A", "B", "H"],
        "genotypes": {
          "A": 0,
          "H": 1,
          "B": 2
        },
        "geno_sep": False,
        "geno_transposed": True,
        "transformations" : []
    }
    write_new_control(control,transformation)

def load_gmap(control):
    """GEMMA2/Rqtl2 eager loading of gmap file"""
    ctrl = methodize(control)
    fn = ctrl.gmap
    logging.info(f"Reading GEMMA2/Rqtl2 gmap {fn}")

def load_geno(control: dict, filter_gs_ok: Callable[[list],bool]) -> np.ndarray:
    """GEMMA2/Rqtl2 eager loading of GENO file. Currently only the compact
format is supported. Returns a genotype numpy matrix with numbers
(minor allele count) and the belonging marker list. filter_gs_ok is a
filtering function for a marker list of genotypes.

    """
    ctrl = methodize(control)
    fn = ctrl.geno
    inds = ctrl.individuals
    markers = ctrl.markers
    genotype_translate = ctrl.genotypes
    print(control)
    assert 'geno_transposed' in control, "Expect geno_transposed set in control file"
    assert markers>inds, f"markers ({markers}) should be larger than individuals ({inds})"
    logging.info(f"Reading GEMMA2/Rqtl2 geno {fn}")
    shape = (markers,inds)
    g = np.empty(shape, dtype=np.float32, order="F")
    # print(g.shape)
    in_header = True
    columns = None
    markerlist = []
    failed = 0
    with gzip.open(fn) as f:
        line = f.readline()
        if line[0] == '#':
            next
        count = None
        while line:
            if in_header:
                in_header = False
                count = 0
                next
            else:
                (marker,l) = line.decode().rstrip().split("\t",3)
                gs = [genotype_translate[v] for v in list(l)]
                assert len(gs)==inds,"number of genotypes for {marker}@{count} differs from {inds}"
                if filter_gs_ok(marker,gs):
                    markerlist.append(marker)
                    g[count,:] = np.array(gs)
                    count += 1
                else:
                    failed += 1
            line = f.readline()
    memory_usage()
    if failed > 0:
        logging.warn(f"{failed} out of {markers} markers were filtered out")
    assert count==markers-failed, f"number of markers ({markers}) does not match ({count+failed}) lines in {fn}"
    return(g,markerlist)

def iter_pheno_txt(fn: str, sep: str = "\t", header: bool = False):
    """Iter of GEMMA2 pheno file. Returns by line"""
    logging.info(f"Reading GEMMA2/Rqtl2 pheno {fn}")
    with open(fn,"r") as f:
        for count, line in enumerate(f):
            if header or count > 1:
                yield count,line.strip().split(sep)

def iter_pheno(fn: str, sep: str = "\t", header: bool = False):
    """Iter of GEMMA2 pheno file. Returns by line"""
    logging.info(f"Reading GEMMA2/Rqtl2 pheno {fn}")
    with gzip.open(fn,"r") as f:
        for count,line in enumerate(f):
            if header or count > 0:
                yield count,line.decode().strip().split(sep)

def iter_geno(fn: str, sep: str = "\t", geno_sep: bool = False, header: bool = False):
    logging.info(f"Reading GEMMA2/Rqtl2 geno {fn}")
    with gzip.open(fn) as f:
        for count,line in enumerate(f):
            if header and count==1:
                h = line.decode()
                hs = h.strip().split(sep)
                yield hs
            if count>1:
                l = line.decode()
                marker,genotypes = l.strip().split(sep,2)
                yield count,marker,[char for char in genotypes]

def iter_geno_num(fn: str, sep: str = "\t", geno_sep: bool = False, header: bool = False):
    for count,marker,gs in iter_geno(fn,sep,geno_sep,header):
        print(marker,gs)
