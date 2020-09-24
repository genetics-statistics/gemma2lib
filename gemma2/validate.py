
import collections
import logging
from os.path import dirname, basename, isfile
import sys

import gemma2.utility.safe as safe
import gemma2.utility.data as data
from gemma2.format.rqtl2 import load_control, iter_pheno, iter_geno, write_new_control

def validate_data(controlfn: dict):
    control = load_control(controlfn)
    ctrl = data.methodize(control)
    path = dirname(controlfn) or "."
    na_strings = ctrl.na_strings
    for num,marker,genos in iter_geno(path+"/"+ctrl.geno, header = False):
        gs = []
        print(genos)
