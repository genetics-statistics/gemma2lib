# Filtering functions

import logging
from os.path import dirname, basename, isfile
import sys
from gemma2.utility.options import get_options_ns
import gemma2.utility.safe as safe
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno, write_control

def filters(controlfn: str, pheno_column: int):
    control = load_control(controlfn)
    print(control)
    na_strings = control.na_strings
    path = dirname(controlfn) or "."
    ids = []
    idx = []
    count = -1

    with safe.pheno_write_open() as p:
        for row in iter_pheno(path+"/"+control.pheno, header = True):
            count += 1
            if not row[pheno_column] in na_strings:
                id = row[0]
                ids.append(id)
                idx.append(count)
                p.write("\t".join([id,row[pheno_column]]).encode())
                p.write("\n".encode())
        phenofn = p.name
    ids = ids[1:] # strip leading column
    idx = idx[1:]
    with safe.geno_write_open() as g:
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
        genofn = g.name
    inds = len(idx)
    # Write control file last
    write_control("filter",inds,control.markers,control.phenotypes,genofn,phenofn)
