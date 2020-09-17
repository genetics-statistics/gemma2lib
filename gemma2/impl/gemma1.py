# GEMMA1 calls

import logging

from gemma2.format.bimbam import write_bimbam
from gemma2.utility.options import get_options_ns

def compute_kinship(control):
    opts = get_options_ns()
    logging.info('Computing GRM with GEMMA1')
    logging.info('Convert to BIMBAM')
    write_bimbam(control.name)
    logging.info('Call gemma')
