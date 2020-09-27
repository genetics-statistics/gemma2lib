# Regression testing
#
# Run with python3 test/regression-tests.py

from glob import glob
from os.path import basename
from pytest_regressions.common import check_text_files
from subprocess import run

run("python3 ./bin/gemma2 --overwrite -o test/data/regression/21487_convert convert --bimbam -g example/21487_BXD_geno.txt.gz  -a example/BXD_snps.txt -p example/21487_BXDPublish_pheno.txt".split())

testdir="test/data/regression"
expectdir=testdir+"/expect"

for efn in glob(expectdir+'/21487_*'):
    if efn.endswith(".gz"):
        continue
    print("Comparing:",efn)
    check_text_files(testdir+"/"+basename(efn),efn)
