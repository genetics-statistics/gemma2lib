# Global options

from types import SimpleNamespace
from gemma2.utility.data import methodize

options = {}

def set_options(opts):
    global options
    options = opts
    options['debug_data'] = options['debug'] == 'DATA' or options['debug'] == 'ALL'
    options['debug_ram'] = options['debug'] == 'RAM' or options['debug'] == 'ALL'

def get_options():
    return options

def get_options():
    return methodize(options)

def get_options_ns():
    return get_options()
