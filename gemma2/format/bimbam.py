# GEMMA2 BIMBAM format support

import json
import gzip
import logging
import numpy as np
from os.path import dirname, basename, splitext, isfile
import sys

from types import SimpleNamespace
from gemma2.utility.options import get_options_ns
from gemma2.utility.system import memory_usage

from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno

# Can not overwrite existing file
class safe_write_open(object):
    def __init__(self, file_name, msg):
        self.file_name = file_name
        self.msg = msg

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: {self.file_name} already exists")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"{self.msg} {self.file_name}")
        self.file.close()

def convert_bimbam(genofn: str, phenofn: str):
    """Read BIMBAM and output to Rqtl2"""
    options = get_options_ns()
    path = options.out_prefix
    outgenofn = path+"_geno.txt.gz"
    logging.info(f"Reading BIMBAM genofile {genofn}")
    logging.info(f"Writing GEMMA2/Rqtl2 genofile {outgenofn}")
    translate = { "1": "A", "0": "B", "0.5": "H" } # FIXME hard coded
    inds = None
    markers = 0
    import gzip
    with gzip.open(outgenofn, mode="wb") as out:
        with gzip.open(genofn, mode='r') as f:
            for line in f:
                markers += 1
                l = line.decode().strip()
                gs = l.split("\t")
                if len(gs) == 1:
                    gs = l.split(", ")
                    assert len(gs) != 1
                out.write(gs[0].encode())
                out.write("\t".encode())
                # print(gs[3:])
                out.write("".join([ translate[item] for item in gs[3:] ]).encode())
                out.write("\n".encode())
                if not inds:
                    inds = len(gs)-3
                    logging.info(f"{inds} individuals")

    logging.info(f"{markers} markers")
    basefn = path
    # Write control file last
    import json
    control = {
        "description": basename(path),
        "crosstype": None,   # we are not assuming a cross for GEMMA
        "sep": "\t",
        "na.strings": ["-", "NA"],
        "comment.char": "#",
        "individuals": inds,
        "markers": markers,
        # "phenotypes": phenos,
        "geno": basename(genofn),
        # "pheno": basename(phenofn),
        "alleles": ["A", "B", "H"],
        "genotypes": {
          "A": 1,
          "H": 2,
          "B": 3
        },
        "geno_sep": False,
        "geno_transposed": True
    }
    controlfn = basefn+".json"
    with safe_write_open(controlfn,"Writing GEMMA2 control file") as f:
        json.dump(control, f, indent=4)

def write_bimbam(controlfn):
    """Write BIMBAM files from R/qtl2 and GEMMA control file"""
    options = get_options_ns()
    path = dirname(controlfn)
    control = load_control(controlfn)
    base = splitext(control.pheno)[0]
    if path:
        base = path + "/" + base

    phenofn = base+"_bimbam.txt"
    logging.info(f"Writing BIMBAM pheno file {phenofn}")
    with open(phenofn,"w") as f:
        for p in iter_pheno(control.pheno, sep=control.sep, header=False):
            # skip the header and the item counter, otherwise same
            f.write("\t".join(p[1:]))

    base = splitext(splitext(control.geno)[0])[0]
    if path:
        base = path + "/" + base
    genofn = base+"_bimbam.txt.gz"
    logging.info(f"Writing BIMBAM geno file {genofn}")
    genotype_translate = { "A": "1", "B": "0", "H": "2"}
    with gzip.open(genofn, mode='wb', compresslevel=options.compression_level) as f:
        # f.write("marker".encode())
        for marker,genotypes in iter_geno(control.geno, sep=control.geno_sep, header=False):
            f.write(marker.encode())
            f.write(" - - ".encode())
            f.write(" ".join([genotype_translate[v] for v in genotypes]).encode())
            f.write("\n".encode())
