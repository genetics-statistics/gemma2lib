import numpy as np
import pandas as pd


def read():
    pass


def read_genotype(filename, mode="bimbam"):
    pass


def read_phenotype(filename):
    pheno = read_tsv_file(filename)
    pheno = pd.DataFrame(pheno)
    return pheno


def read_covariates(filename):
    covariates = read_tsv_file(filename)
    covariates = pd.DataFrame(covariates)
    return covariates


def read_control():
    pass


def read_tsv_file(filename):
    data = []
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break
            data.append(str(line).rstrip().split("\t"))
    return data
