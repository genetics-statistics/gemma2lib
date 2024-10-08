# GEMMA2 BIMBAM format support

import copy
import json
import gzip
import logging
import numpy as np
from os.path import dirname, basename, splitext, isfile
import re
import sys

import gemma2.utility.safe as safe
from types import SimpleNamespace
from gemma2.utility.options import get_options_ns
from gemma2.utility.system import memory_usage

from gemma2.format.rqtl2 import load_control, write_control, iter_pheno, iter_geno, methodize

def convert_bimbam(genofn: str, phenofn: str, annofn: str):
    """Read/convert/import BIMBAM and output to Rqtl2 format"""
    options = get_options_ns()
    path = options.out_prefix

    basefn = path

    logging.info(f"Reading BIMBAM marker/SNP {annofn}")
    with open(annofn,"r") as f:
        with safe.gmap_write_open() as out:
            outgmapfn = out.name
            out.write(f"marker,chr,pos\n".encode())
            for line in f:
                list = re.split('[,\t\s]+',line.strip())
                # print(list)
                marker,pos,chr = list
                out.write(f"{marker}\t{chr}\t{pos}\n".encode())

    logging.info(f"Reading BIMBAM phenofile {phenofn}")
    in_header = True
    p_inds = 0
    phenos = None
    with open(phenofn,"r") as f:
        with safe.pheno_write_open("_pheno_bimbam.txt.gz") as out:
            outphenofn = out.name
            for line in f:
                ps = line.strip().split("\t")
                if not phenos:
                    phenos = len(ps)
                if in_header:
                    out.write("id\t".encode())
                    out.write("\t".join([f"{i+1}" for i in range(phenos)]).encode())
                    in_header = False
                    out.write("\n".encode())
                p_inds += 1
                out.write(f"{p_inds}\t".encode())
                out.write("\t".join(ps).encode())
                out.write("\n".encode())

    inds = None
    markers = 0
    translate = { "1": "A", "0": "B", "0.5": "H" } # FIXME hard coded

    in_header = True
    with safe.geno_write_open("_geno_bimbam.txt.gz") as out:
        outgenofn = out.name
        with gzip.open(genofn, mode='r') as f:
            for line in f:
                markers += 1
                l = line.decode().strip()
                gs = l.split("\t")
                if len(gs) == 1:
                    gs = l.split(", ")
                genos = len(gs)-3
                assert genos != 1
                if in_header:
                    out.write("marker\t".encode())
                    out.write("\t".join([f"{i+1}" for i in range(genos)]).encode())
                    in_header = False
                    out.write("\n".encode())
                out.write(gs[0].encode())
                out.write("\t".encode())
                # print(gs[3:])
                out.write("".join([ translate[item] for item in gs[3:] ]).encode())
                out.write("\n".encode())
                if not inds:
                    inds = len(gs)-3
                    logging.info(f"{inds} individuals")

    logging.info(f"{markers} markers")
    logging.info(f"{phenos} phenotypes")
    assert inds == p_inds, f"Individuals not matching {inds} != {p_inds}"
    transformation = { "type": "export", "original": "rqtl2", "format": "bimbam" }
    write_control(None,inds,markers,phenos,outgenofn,outphenofn,outgmapfn,transformation)

def write_bimbam(control: dict):
    """Write BIMBAM files from R/qtl2 and GEMMA control file"""
    ctrl = methodize(control)

    with safe.pheno_write_open("_pheno_bimbam.txt") as f:
        phenofn = f.name
        for num,p in iter_pheno(ctrl.pheno, sep=ctrl.sep, header=False):
            # skip the header and the item counter, otherwise same
            f.write("\t".join(p[1:]))
            f.write("\n")

    genotype_translate = ctrl.genotypes
    # This provides swapped values for major allele swap
    genotype_translate_swap = copy.copy(ctrl.genotypes)
    a = copy.copy(genotype_translate['A'])
    b = copy.copy(genotype_translate['B'])
    genotype_translate_swap['A'] = b
    genotype_translate_swap['B'] = a

    genoA = ctrl.alleles[0]
    genoB = ctrl.alleles[1]
    with safe.geno_write_open("_bimbam.txt.gz") as f:
        genofn = f.name
        for num,marker,genotypes in iter_geno(ctrl.geno, sep=ctrl.sep, geno_sep=ctrl.geno_sep, header=False):
            f.write(marker.encode())
            # we need to write major alleles - GEMMA expects them to be 0.0
            histogram = { 'A': 0, 'B': 0 }
            for v in genotypes:
                if v not in histogram:
                    histogram[v] = 0
                histogram[v] += 1
            # Is A the real major allele?
            if histogram['A'] > histogram['B']:
                translate = genotype_translate
                f.write(f",{genoB},{genoA},".encode()) # minor allele first
            else:
                translate = genotype_translate_swap
                f.write(f",{genoA},{genoB},".encode()) # minor allele first

            values = [str(translate[v]) if v!="-" else "NA" for v in genotypes]
            f.write(",".join(values).encode())
            f.write("\n".encode())
    return genofn, phenofn
