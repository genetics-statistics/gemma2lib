import logging
from os.path import dirname, basename
from types import SimpleNamespace

from pandas_plink import read_plink

import os, psutil, numpy as np
def memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_full_info()[0]
    if (mem > 1024**3):
        print(f"{round(mem / float(2 ** 20) / 102.4 )/10}Gb RAM used")
    else:
        print(f"{round(mem / float(2 ** 20) )}Mb RAM used")

def convert_plink(path: str, options: SimpleNamespace):
    """Convert PLINK format to GEMMA2"""
    verbose = options.verbose
    (bim,fam,bed) = read_plink(path, verbose=(True if verbose>1 else False))
    m = bed.compute()
    if options.debug:
        print("Debug view of PLINK\n")
        print("===> BIM") # alleles/markers
        print(bim.head())
        print(bim.info())
        print("===> FAM") # samples/phenotypes
        print(fam.head())
        print(fam.info())
        print("===> BED") # genotypes
        print(m)
        print(m[0:3,0:40])

    print(bim.shape)
    markers2,phenos2 = bim.shape
    inds2,phenos = fam.shape
    markers, inds = m.shape
    assert inds == inds2, "Number of individuals not matching in fam and bed files"
    assert markers == markers2, "Number of markers not matching in bim and bed files"
    assert phenos == phenos2, "Number of phenotypes not matching in bim and fam files"
    if options.debug or options.verbose>2:
        memory_usage()
    # Writing genotype file
    genofn = basename(path)+"_geno.tsv"
    logging.info(f"Writing geno file {genofn}")
    translate = { 1.0: "A", 2.0: "B", 0.0: "H" }
    with open(genofn, "w") as f:
        f.write("marker")
        for i in range(inds):
            f.write(f"\t{i}")
        for j in range(markers):
            markername = bim.snp[j]
            f.write(f"\n{markername}")
            for i in range(inds):
                f.write(f"\t{translate[m[j,i]]}")
