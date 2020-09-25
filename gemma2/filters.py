# Filtering functions

import collections
import logging
from os.path import dirname, basename, isfile
import sys
from gemma2.utility.options import get_options_ns
import gemma2.utility.safe as safe
import gemma2.utility.data as data
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno, write_new_control

def error(msg: str):
    raise ValueError(msg)

def maf_filter(marker: str, maf_threshold: float, miss_threshold: float,
               gs: list, na_strings: list) -> bool:
    """Pass maf threshold? FIXME: need to account for values.
    Returns True if SNP passes MAF threshold
    """
    num = len(gs)
    realgs = list(filter(lambda g: not g in na_strings, gs))
    realnum = len(realgs)
    missing = num - realnum
    counter=collections.Counter(realgs)
    miss_fract = missing/num
    if miss_fract > miss_threshold:
        logging.debug(f"Miss filter {miss_fract} drops {counter} at {marker}")
        return False
    if len(counter) < 2:
        logging.debug(f"MAF filter {maf_threshold} drops {counter} at {marker}")
        return False
    # we take the second value which differs from GEMMA1 in the rare
    # instance that we have enough Heterozygous - FIXME when we have
    # genotype numbers - H should count by minor allele 50%.
    minor_count = counter.most_common()[1][1]
    ret = minor_count/realnum > maf_threshold
    if not ret:
        logging.debug(f"MAF filter {maf_threshold} drops {counter} at {marker}")
    return ret

def filters(controlfn: str, pheno_column: int, maf: float, miss: float):
    control = load_control(controlfn)
    ctrl = data.methodize(control)
    na_strings = ctrl.na_strings
    path = dirname(controlfn) or "."
    ids = []
    idx = []
    count = -1 # FIXME, now comes with iter

    with safe.pheno_write_open() as p:
        logging.debug(f"Filtering phenotypes on missing data (NA)")
        for num,row in iter_pheno(path+"/"+ctrl.pheno, header = True):
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
    markers = 0
    with safe.geno_write_open() as g:
        logging.debug(f"Filtering genotypes on missing data at {miss}")
        logging.debug(f"Filtering genotypes on maf at {maf}")
        g.write("\t".join(["marker"]+ids).encode())
        g.write("\n".encode())
        for num,marker,genos in iter_geno(path+"/"+ctrl.geno, header = False):
            gs = []
            for i in idx:
                gs.append(genos[i-1])
            if maf_filter(marker,maf,miss,gs,na_strings):
                markers += 1
                g.write(marker.encode())
                g.write("\t".encode())
                g.write("".join(gs).encode())
                g.write("\n".encode())
        genofn = g.name
    inds = len(idx)
    # Write control file last
    import copy
    ncontrol = copy.deepcopy(control)
    ncontrol['command'] = "filter"
    ncontrol['geno'] = genofn
    ncontrol['geno_compact'] = True
    ncontrol['geno_transposed'] = True
    ncontrol['pheno'] = phenofn
    ncontrol['individuals'] = inds
    ncontrol['markers'] = markers
    ncontrol['maf'] = maf
    write_new_control(ncontrol)
