import os, psutil, numpy as np

def memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_full_info()[0]
    if (mem > 1024**3):
        print(f"{round(mem / float(2 ** 20) / 102.4 )/10}Gb RAM used")
    else:
        print(f"{round(mem / float(2 ** 20) )}Mb RAM used")
