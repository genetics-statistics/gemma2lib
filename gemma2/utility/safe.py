import logging
from os.path import dirname, basename, isfile
import sys
from gemma2.utility.options import get_options_ns

# Can not overwrite existing file
import gzip
import json

class control_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix+".json"

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: control file {self.file_name} already exists")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 control {self.file_name}")
        self.file.close()

class pheno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix+"_pheno.txt.gz"
        self.compression_level = opts.compression_level

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'wb', compresslevel=self.compression_level)
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 geno {self.file_name}")
        self.file.close()

class geno_write_open(object):
    def __init__(self):
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix+"_geno.txt.gz"
        self.compression_level = opts.compression_level

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'wb', compresslevel=self.compression_level)
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 geno {self.file_name}")
        self.file.close()
