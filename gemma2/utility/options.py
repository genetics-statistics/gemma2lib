# Global options

from types import SimpleNamespace

options = {}

def set_options(opts):
    global options
    options = opts
    options['debug_data'] = options['debug'] == 'DATA' or options['debug'] == 'ALL'
    options['debug_ram'] = options['debug'] == 'RAM' or options['debug'] == 'ALL'

def get_options():
    return options

def get_options_ns():
    return SimpleNamespace(**options)
