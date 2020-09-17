# GRM/Kinship

from gemma2.format.rqtl2 import load_control, load_geno
from gemma2.utility.options import get_options_ns
import gemma2.impl.gemma1 as gemma1

def compute_kinship(control,impl,loco):
    """GRM/Kinship computation"""
    opts = get_options_ns()
    if impl == "gemma1":
        print(control)
        gemma1.compute_kinship(control)
    else:
        print(control)
        G = load_geno(control)
        print(G)
