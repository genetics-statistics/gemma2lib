# Global options

from types import SimpleNamespace

options = {}

def set_options(opts):
    global options
    options = opts

def get_options():
    return options

def get_options_ns():
    return SimpleNamespace(**options)
