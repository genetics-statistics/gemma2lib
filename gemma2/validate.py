
import collections
import logging
from os.path import dirname, basename, isfile
import sys

import gemma2.utility.safe as safe
import gemma2.utility.data as data
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

def warn_once(id: str, msg: str):
    global _warnings
    if not id in _warnings:
        _warnings[id] = True
        msg += f" --- other similar {id} warnings are ignored"
        warn(msg)

def error(msg: str):
    msg += f" ({_err_msg} file {_err_filename} line {_err_line})"
    logging.error(msg)
    raise ValidationError(msg)

def check_genotypes(num,gs,na_strings):
    # print(gs)
    realgs = list(filter(lambda g: not g in na_strings, gs))
    realnum = len(realgs)
    fract_real = realnum/len(gs)
    if fract_real < CONST_FRACT_REAL:
        warn(f"Low fraction of real genotypes ({fract_real}) in {gs}")
    if realnum < CONST_MIN_REAL:
        warn_once("realnum",f"There are not enough genotypes ({realnum}) in {gs}")
    missing = len(gs) - realnum
    counter=collections.Counter(realgs)
    if len(counter) < 2:
        warn_once("counter",f"Minimal two genotypes required {counter} found in {gs}")
    # we take the second value which differs from GEMMA1 in the rare
    # instance that we have enough Heterozygous - FIXME when we have
    # genotype numbers - H should count by minor allele 50%.
    minor_count = counter.most_common()[1][1]
    # ret = minor_count/realnum > maf_threshold
    # if not ret:
    #     logging.debug(f"MAF filter {maf_threshold} fails {counter} at {marker}")
    # return ret


def validate_data(controlfn: dict):
    control = load_control(controlfn)
    ctrl = data.methodize(control)
    path = dirname(controlfn) or "."
    na_strings = ctrl.na_strings
    fn = path + "/" + ctrl.geno
    global _err_msg, _err_filename, _err_line
    _err_filename = fn
    for num,marker,genos in iter_geno(fn, header = False):
        _err_msg = marker
        _err_line = num
        check_genotypes(num,genos,na_strings)
