import os, psutil, numpy as np
from gemma2.utility.options import get_options_ns

def memory_usage(msg: str = None):
    options = get_options_ns()
    if options.debug_ram or options.verbose>2:
        process = psutil.Process(os.getpid())
        mem = process.memory_full_info()[0]
        if msg:
            print(msg,end=":\t")
        if (mem > 1024**3):
            print(f"{round(mem / float(2 ** 20) / 102.4 )/10}Gb RAM used")
        else:
            print(f"{round(mem / float(2 ** 20) )}Mb RAM used")
