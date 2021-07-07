import pytest
import pandas as pd

from gemma2.io import read_phenotype, read_covariates


def test_read_phenotype():
    phenotype = read_phenotype("./example/mouse_hs1940.pheno.txt")
    assert isinstance(phenotype, pd.DataFrame)
    assert phenotype.shape == (1940, 6)


def test_read_covariates():
    covariates = read_covariates("./example/BXD_covariates.txt")
    assert isinstance(covariates, pd.DataFrame)
    assert covariates.shape == (198, 2)
