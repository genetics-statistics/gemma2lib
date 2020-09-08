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
        logging.info(f"Writing GEMMA2/Rqtl2 pheno {self.file_name}")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        self.file.close()

# Can not overwrite existing file
class safe_geno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.file_name = opts.out_prefix+"_geno.txt"

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        logging.info(f"Writing GEMMA2/Rqtl2 geno {self.file_name}")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
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
        g.write("\t".join(["marker"]+ids[1:]))
        g.write("\n")
        for marker,genos in iter_geno(path+"/"+control.geno, header = False):
            gs = []
            for i in idx:
                gs.append(genos[i-1])
            g.write(marker)
            g.write("\t")
            g.write("".join(gs))
            g.write("\n")
