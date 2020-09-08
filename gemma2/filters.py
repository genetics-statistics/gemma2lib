# Filtering functions

import logging
from os.path import dirname, basename, isfile
from gemma2.utility.options import get_options_ns
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno

# Can not overwrite existing file
class safe_pheno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.file_name = opts.out_prefix+"_pheno.txt"

    def __enter__(self):
        if isfile(self.file_name):
            raise Exception(f"ERROR: {self.file_name} already exists")
        logging.info(f"Open pheno file {self.file_name}")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        self.file.close()

def filter(controlfn: str, pheno_column: int):
    control = load_control(controlfn)
    print(control)
    na_strings = control.na_strings
    path = dirname(controlfn)
    with safe_pheno_write_open() as p:
        for row in iter_pheno(path+"/"+control.pheno, header = True):
            if not row[pheno_column] in na_strings:
                p.write("\t".join([row[0],row[pheno_column]]))
                p.write("\n")
