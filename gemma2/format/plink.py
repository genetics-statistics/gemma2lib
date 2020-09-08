import logging
from os.path import dirname, basename
import sys

from gemma2.utility.options import get_options_ns
from pandas_plink import read_plink
from gemma2.utility.system import memory_usage

def convert_plink(path: str, compression_level: int):
    """Convert PLINK format to GEMMA2"""
    def mknum(v):
        if v != v:
            return "NA"
        return(str(v))

    options = get_options_ns()
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

    phenofn = basefn+"_pheno.tsv"
    logging.info(f"Writing GEMMA2 pheno file {phenofn}")
    p = fam.to_numpy()
    with open(phenofn, mode="w") as f:
        f.write("id")
        for c in fam.columns.values:
            if c != "i": # we skip the last i column
                f.write(f"\t{c}")
        f.write("\n")
        for j in range(inds):
            f.write(str(j+1)+"\t")
            f.write("\t".join([mknum(v) for v in p[j,:-1]])) # except for i column
            f.write("\n")

    memory_usage("plink pheno")

    genofn = basefn+"_geno.tsv.gz"
    logging.info(f"Writing GEMMA2 geno file {genofn}")
    translate = { 1.0: "A", 2.0: "B", 0.0: "H" }

    import gzip
    # content = b"Lots of content here"
    with gzip.open(genofn, mode='wb', compresslevel=compression_level) as f:
    # with open(genofn, "w") as f:
        f.write("marker".encode())
        for i in range(inds):
            f.write(f"\t{i}".encode())
        for j in range(markers):
            markername = bim.snp[j]
            f.write(f"\n{markername}\t".encode())
            if options.low_mem: # shaves 20%
                for i in range(inds):
                    f.write(f"{translate[m[j,i]]}".encode())
            else:
                f.write("".join([ translate[item] for item in m[j] ]).encode())

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
        "phenotypes": phenos,
        "geno": basename(genofn),
        "pheno": basename(phenofn),
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
    logging.info(f"Writing GEMMA2 control file {controlfn}")
    with open(controlfn, 'w') as cf:
        json.dump(control, cf, indent=4)

    memory_usage("plink geno")
