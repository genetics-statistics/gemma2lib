import logging
from os.path import dirname, basename, isfile, exists
import sys
import numpy as np

from gemma2.utility.options import get_options_ns
import gemma2.utility.safe as safe
from pandas_plink import read_plink
from gemma2.utility.system import memory_usage

from gemma2.format.rqtl2 import write_control

def convert_plink(path: str, annofn: str):
    """Convert PLINK format to GEMMA2"""
    def mknum(v):
        if v != v:
            return "NA"
        return(str(v))

    options = get_options_ns()
    compression_level = options.compression_level
    verbose = options.verbose
    memory_usage("plink before load")

    logging.info(f"Reading PLINK files {path}")
    (bim,fam,bed) = read_plink(path, verbose=(True if verbose>1 else False))
    m = bed.compute()
    if options.debug_data:
        print("Debug view of PLINK\n")
        print("===> BIM alleles/markers")
        print(bim.head())
        print(bim.info())
        print("===> FAM samples/phenotypes")
        print(fam.head())
        print(fam.info())
        print("===> BED genotypes")
        print(m)
        print([x for x in m[0]])
        print(bim.shape)

    markers2,phenos2 = bim.shape
    inds2,phenos = fam.shape
    markers, inds = m.shape
    assert inds == inds2, "Number of individuals not matching in fam and bed files"
    assert markers == markers2, "Number of markers not matching in bim and bed files"
    assert phenos == phenos2, "Number of phenotypes not matching in bim and fam files"

    basefn = options.out_prefix
    memory_usage("plink pandas")

    outgmapfn = "NONE"
    if annofn:
        logging.info(f"Reading BIMBAM marker/SNP {annofn}")
        with open(annofn,"r") as f:
            with safe.gmap_write_open() as out:
                outgmapfn = out.name
                out.write(f"marker,chr,pos\n".encode())
                for line in f:
                    marker,pos,chr,rest = line.strip().split("\t")
                    out.write(f"{marker}\t{chr}\t{pos}\n".encode())

    phenofn = basefn+"_pheno.tsv"
    p = fam.to_numpy()
    with safe.pheno_write_open() as f:
        outphenofn = f.name
        f.write("id".encode())
        for c in fam.columns.values:
            if c != "i": # we skip the last i column
                f.write(f"\t{c}".encode())
        f.write("\n".encode())
        for j in range(inds):
            f.write((str(j+1)+"\t").encode())
            f.write("\t".join([mknum(v) for v in p[j,:-1]]).encode()) # except for i column
            f.write("\n".encode())

    memory_usage("plink pheno")

    genofn = basefn+"_geno.txt.gz"
    logging.info(f"Writing GEMMA2 geno file {genofn}")
    translate = { 1.0: "A", 2.0: "B", 0.0: "H", -9.0: "-" }

    import gzip
    with gzip.open(genofn, mode='wb', compresslevel=compression_level) as f:
        f.write("marker".encode())
        for i in range(inds):
            f.write(f"\t{i+1}".encode())
        for j in range(markers):
            markername = bim.snp[j]
            f.write(f"\n{markername}\t".encode())
            values = [ -9.0 if np.isnan(x) else x for x in m[j] ]
            if options.low_mem: # shaves 20%
                for i in range(inds):
                    f.write(f"{translate[values[i]]}".encode())
            else:
                f.write("".join([ translate[item] for item in values ]).encode())
        outgenofn = genofn

    transformation = { "type": "convert", "original": "plink", "format": "rqtl2" }
    write_control(None,inds,markers,phenos,outgenofn,outphenofn,outgmapfn,transformation)

    memory_usage("plink geno")
