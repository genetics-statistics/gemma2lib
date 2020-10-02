# Run doctest with
# env PYTHONPATH=$PYTHONPATH:. python3 -m doctest -v gemma2/compute/geno.py

import collections
import logging

def real_genotypes(gs: list, na_strings: dict) -> list:
    """Strips out None and NAs"""
    return list(filter(lambda g: g!=None and not g in na_strings, gs))

def real_alleles(gs: list, na_strings: dict, h_str:str) -> list:
    """Returns major and minor alleles only"""
    return list(filter(lambda g: g!=None and g!=h_str and not g in na_strings, gs))

def calc_miss(gs: list, na_strings: dict) -> float:
    num = len(gs)
    realgs = real_genotypes(gs,na_strings)
    realnum = len(realgs)
    missing = num - realnum
    counter=collections.Counter(realgs)
    return missing/num

def is_miss_fail(gs: list, miss_threshold: float, na_strings: dict) -> bool:
    return calc_miss(gs,na_strings) > miss_threshold

def calc_maf(gs: list, na_strings: dict, h_str: str) -> float:
    """
    Compute minor allele frequency.

    >>> calc_maf(["A","A","A","B","B","B","A","H"],{"NA": 1},"H")
    0.4375

    >>> calc_maf(["B","B","B","A","A","A","B","H"],{"NA": 1},"H")
    0.4375

    >>> calc_maf(["B","B","B","A","A","A","B","H","NA"],{"NA": 1},"H")
    0.4375

    >>> calc_maf(["B","B","B","A","A","A","B","H","H","H","NA"],{"NA": 1},"H")
    0.45

    >>> calc_maf(["B","B","B","A","A","A","B","H","H","H","H","NA"],{"NA": 1},"H")
    0.45454545454545453

    """
    values = list(filter(lambda g: g!=None,geno_translate_to_num(gs,na_strings,h_str)))
    # print(values)
    realnum = len(values)
    # print(len(values),values)
    minor_count = 0.0
    for geno in values:
        minor_count += geno/2.0
    # print(minor_count,realnum)
    return minor_count/realnum

def is_maf_fail(gs: list, maf_threshold: float, na_strings: dict) -> bool:
    return calc_maf(gs,na_strings) < maf_threshold

def geno_translate_to_num(gs: list, na_strings: dict = {"NA": 1}, h_str: str = "H") -> list:
    """Translate (3 char) genotypes to minor allele frequencies

    >>> geno_translate_to_num(["A","A","A","B","B","H"])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0]

    >>> geno_translate_to_num(["A","A","A","B","B","H",None])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0, None]

    Or with NA

    >>> geno_translate_to_num(["A","A","A","B","B","H","NA"])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0, None]

    We can reverse alleles and get the minor count

    >>> geno_translate_to_num(["B","B","B","A","A","H",None])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0, None]

    >>> geno_translate_to_num(["B","B","B","A","H","H",None])
    [0.0, 0.0, 0.0, 2.0, 1.0, 1.0, None]

    >>> geno_translate_to_num(["B","B","B","A","H","H",None,"NA"])
    [0.0, 0.0, 0.0, 2.0, 1.0, 1.0, None, None]

    These should fail

    >>> geno_translate_to_num(["B","B","B","B",None]) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    AssertionError: Expected multiple alleles for Counter({'B': 4})

    >>> geno_translate_to_num(["B","B","B","B","A","H","C"]) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    AssertionError: Expected max 3 alleles for Counter({'B': 4, 'A': 1, 'H': 1, 'C': 1})

    """
    realgs = real_alleles(gs,na_strings,h_str)
    counter=collections.Counter(realgs)
    count = counter.most_common()
    alleles = len(count)
    assert alleles>1, f"Expected multiple alleles for {counter}"
    assert alleles<=2, f"Expected max 2 alleles for {counter}"
    trans = {}
    trans[count[0][0]] = 0.0 # Major allele
    if alleles > 1:
        trans[count[1][0]] = 2.0 # Set minor allele
    trans[h_str] = 1.0 # Set heterozygous
    trans[None] = None
    for k in na_strings:
        trans[k] = None
    return [ trans[g] for g in gs ]
