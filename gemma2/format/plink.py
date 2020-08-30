import logging

from pandas_plink import read_plink

def convert_plink(path: str, verbose: int, debug: bool):
    (bim,fam,bed) = read_plink(path, verbose=(True if verbose>1 else False))
    if debug:
        print("Debug view of PLINK\n")
        print("===> BED") # alleles/markers
        print(bim.head())
        print("===> BIM") # samples/phenotypes
        print(fam.head())
        print("===> BED") # genotypes
        m = bed.compute()
        print(m)
        print(m[0:3,0:40])
