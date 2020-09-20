import logging
from os.path import dirname, basename, isfile
import sys
from gemma2.utility.options import get_options_ns

# Can not overwrite existing file
import gzip
import json

class write_open:
    def __init__(self, file_type: str, postfix: str = None):
        """file_type is control"""
        self.file_type = file_type
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix
        if postfix:
            self.file_name += postfix
        self.compression_level = opts.compression_level

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: {self.file_type} file {self.file_name} already exists")
        self.file = open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 {self.file_type} {self.file_name}")
        self.file.close()

class control_write_open(write_open):
    def __init__(self):
        write_open.__init__(self,"control",".json")

class pheno_write_open:
    def __init__(self):
        write_open.__init__(self,"pheno","_pheno.txt.gz")

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: pheno file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'wb', compresslevel=self.compression_level)
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 {self.file_type} {self.file_name}")
        self.file.close()

class geno_write_open:
    def __init__(self, postfix: str = "_geno.txt.gz"):
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix+postfix
        self.compression_level = opts.compression_level

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: geno file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'wb', compresslevel=self.compression_level)
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 geno {self.file_name}")
        self.file.close()

class gmap_write_open:
    def __init__(self):
        opts = get_options_ns()
        self.opts = opts
        self.file_name = opts.out_prefix+"_gmap.txt.gz"
        self.compression_level = opts.compression_level

    def __enter__(self):
        if not self.opts.overwrite and isfile(self.file_name):
            raise Exception(f"ERROR: marker/SNP file {self.file_name} already exists")
        self.file = gzip.open(self.file_name, 'w')
        return self.file

    def __exit__(self, type, value, tb):
        logging.info(f"Writing GEMMA2/Rqtl2 marker/SNP {self.file_name}")
        self.file.close()
