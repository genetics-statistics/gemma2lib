# Filtering functions

import logging
from os.path import dirname, basename, isfile
import sys
from gemma2.utility.options import get_options_ns
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno

# Can not overwrite existing file
class safe_pheno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.file_name = opts.out_prefix+"_pheno.txt"

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 pheno {self.file_name}")
        self.file.close()

# Can not overwrite existing file
import gzip
class safe_geno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.file_name = opts.out_prefix+"_geno.txt.gz"
        self.compression_level = opts.compression_level

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'wb', compresslevel=self.compression_level)
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 geno {self.file_name}")
        self.file.close()

import json

class safe_control_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.file_name = opts.out_prefix+".json"

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 control {self.file_name}")
        self.file.close()


def filter(controlfn: str, pheno_column: int):
    control = load_control(controlfn)
    print(control)
    na_strings = control.na_strings
    path = dirname(controlfn) or "."
    ids = []
    idx = []
    count = -1

    with safe_pheno_write_open() as p:
        for row in iter_pheno(path+"/"+control.pheno, header = True):
            count += 1
            if not row[pheno_column] in na_strings:
                id = row[0]
                ids.append(id)
                idx.append(count)
                p.write("\t".join([id,row[pheno_column]]))
                p.write("\n")
    ids = ids[1:] # strip leading column
    idx = idx[1:]
    with safe_geno_write_open() as g:
        g.write("\t".join(["marker"]+ids).encode())
        g.write("\n".encode())
        for marker,genos in iter_geno(path+"/"+control.geno, header = False):
            gs = []
            for i in idx:
                gs.append(genos[i-1])
            g.write(marker.encode())
            g.write("\t".encode())
            g.write("".join(gs).encode())
            g.write("\n".encode())
    inds = len(idx)
    # Write control file last
    opts = get_options_ns()
    control = {
        "description": basename(path),
        "crosstype": None,   # we are not assuming a cross for GEMMA
        "sep": "\t",
        "na.strings": ["-"],
        "comment.char": "#",
        "individuals": inds,
        "markers": control.markers,
        "phenotypes": control.phenotypes,
        "geno": opts.out_prefix+"_geno.txt.gz",
        "pheno": opts.out_prefix+"_pheno.txt",
        "alleles": ["A", "B", "H"],
        "genotypes": {
          "A": 0,
          "H": 1,
          "B": 2
        },
        "geno_sep": False,
        "geno_transposed": True
    }
    with safe_control_write_open() as controlf:
        json.dump(control, controlf, indent=4)
