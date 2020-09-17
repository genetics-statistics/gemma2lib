# GRM/Kinship

from gemma2.format.rqtl2 import load_control, load_geno

def compute_kinship(control,impl,loco):
    """GRM/Kinship computation"""
    print(control)
    G = load_geno(control)
    print(G)
