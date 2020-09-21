# GRM/Kinship

from gemma2.format.rqtl2 import load_control, load_geno
from gemma2.utility.options import get_options_ns
import gemma2.impl.gemma1 as gemma1grm
import gemma2.impl.gemma2.kinship as gemma2grm

def compute_kinship(control,impl,loco):
    """GRM/Kinship computation"""
    opts = get_options_ns()
    if impl == "gemma1":
        gemma1grm.compute_kinship(control)
    else:
        gemma2grm.compute_kinship(control)
