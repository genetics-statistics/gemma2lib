
import collections
import logging
from os.path import dirname, basename, isfile
import sys

import gemma2.utility.safe as safe
import gemma2.utility.data as data
from gemma2.compute.geno import is_maf_fail, is_miss_fail
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno, write_new_control

CONST_FRACT_REAL = 0.3  # allow 70% missing genotype data on a marker
CONST_MIN_REAL   = 10   # minimum 10 real genotypes on a marker

# Using module vars to keep track
_err_msg=""
_err_filename=None
_err_line=-1
_warnings = {}

class ValidationError(ValueError):
    pass

def warn(msg: str):
    msg += f" ({_err_msg} file {_err_filename} line {_err_line})"
    logging.warn(msg)

def warn_limit(id: str, msg: str):
    """Limit the number of warnings"""
    global _warnings
    if not id in _warnings:
        _warnings[id] = 0
    _warnings[id] += 1
    if _warnings[id] > 3:
        return
    if _warnings[id] == 3:
        msg += f" --- other similar {id} warnings are ignored"
    warn(msg)

def error(msg: str):
    msg += f" ({_err_msg} file {_err_filename} line {_err_line})"
    logging.error(msg)
    raise ValidationError(msg)

def check_genotypes(num: int, gs: list, na_strings: dict, maf: float,miss: float) -> dict:
    realgs = list(filter(lambda g: not g in na_strings, gs))
    realnum = len(realgs)
    fract_real = realnum/len(gs)
    if fract_real < CONST_FRACT_REAL:
        warn(f"Low fraction of real genotypes ({fract_real}) in {gs}")
    if realnum < CONST_MIN_REAL:
        warn_limit("realnum",f"There are not enough genotypes ({realnum}) in {gs}")
    missing = len(gs) - realnum
    counter=collections.Counter(realgs)
    if len(counter) < 2:
        warn_limit("counter",f"Minimal two genotypes required {counter} found in {gs}")

    if is_miss_fail(gs,miss,na_strings):
        counter["miss_fail"] = 1
    if is_maf_fail(gs,maf,na_strings):
        counter["maf_fail"] = 1
    return counter

def validate_data(controlfn: dict, maf: float, miss: float):
    control = load_control(controlfn)
    ctrl = data.methodize(control)
    path = dirname(controlfn) or "."
    na_strings = ctrl.na_strings
    fn = path + "/" + ctrl.geno
    global _err_msg, _err_filename, _err_line
    _err_filename = fn
    histogram = {}
    marker_allele_count = {}
    for num,marker,genos in iter_geno(fn, header = False):
        _err_msg = marker
        _err_line = num
        counts = check_genotypes(num,genos,na_strings,maf,miss)
        for item in counts:
            if item not in histogram:
                histogram[item] = 0
            histogram[item] += counts[item]
        # ---- Track counts of major allele by marker:
        major_allele = counts.most_common()[0][0]
        if major_allele not in marker_allele_count:
            marker_allele_count[major_allele] = 0
        marker_allele_count[major_allele] += 1
    markers = num-1
    inds = len(genos)

    if marker_allele_count[major_allele] == markers:
        warn(f"The same major allele is used for all markers; alleles {marker_allele_count}")

    logging.info(f"In {fn} counted {inds} markers {markers} samples and {histogram}")
    assert ctrl.markers == markers
    assert ctrl.individuals == inds
    if len(_warnings) == 0:
        logging.info(f"Great! Genotype file {fn} looks OK")
    else:
        logging.warn(f"Genotype file {fn} has warnings")
        sys.exit(2)
