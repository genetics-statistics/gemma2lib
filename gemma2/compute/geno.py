# Run doctest with
# env PYTHONPATH=$PYTHONPATH:. python3 -m doctest -v gemma2/compute/geno.py

import collections
import logging

def calc_miss(gs: list, na_strings: dict) -> float:
    num = len(gs)
    realgs = list(filter(lambda g: not g in na_strings, gs))
    realnum = len(realgs)
    missing = num - realnum
    counter=collections.Counter(realgs)
    return missing/num

def is_miss_fail(gs: list, miss_threshold: float, na_strings: dict) -> bool:
    return calc_miss(gs,na_strings) > miss_threshold

def calc_maf(gs: list, na_strings: dict) -> float:
    """
    >>> calc_maf(["A","A","A","B","B","B","A","H"],{"NA": 1})
    0.375

    >>> calc_maf(["B","B","B","A","A","A","B","H"],{"NA": 1})
    0.375

    >>> calc_maf(["B","B","B","A","A","A","B","H","NA"],{"NA": 1})
    0.375

    """
    num = len(gs)
    realgs = list(filter(lambda g: not g in na_strings, gs))
    realnum = len(realgs)
    missing = num - realnum
    counter=collections.Counter(realgs)
    # we take the second value which differs from GEMMA1 in the rare
    # instance that we have enough Heterozygous - FIXME when we have
    # genotype numbers - H should count by minor allele 50%.
    minor_count = counter.most_common()[1][1]
    return minor_count/realnum

def is_maf_fail(gs: list, maf_threshold: float, na_strings: dict) -> bool:
    return calc_maf(gs,na_strings) < maf_threshold

def geno_translate_to_num(gs: list) -> list:
    """Translate (3 char) genotypes to minor allele frequencies

    >>> geno_translate_to_num(["A","A","A","B","B","H"])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0]

    >>> geno_translate_to_num(["A","A","A","B","B","H",None])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0, None]

    We can reverse alleles and get the minor count

    >>> geno_translate_to_num(["B","B","B","A","A","H",None])
    [0.0, 0.0, 0.0, 2.0, 2.0, 1.0, None]

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
    values = filter(None, gs)
    counter=collections.Counter(values)
    # we take the second value which differs from GEMMA1 in the rare
    # instance that we have enough Heterozygous - FIXME when we have
    # genotype numbers - H should count by minor allele 50%.
    # print(counter)
    count = counter.most_common()
    alleles = len(count)
    assert alleles>1, f"Expected multiple alleles for {counter}"
    assert alleles<=3, f"Expected max 3 alleles for {counter}"
    trans = {}
    trans[count[0][0]] = 0.0 # Major allele
    if alleles > 1:
        trans[count[1][0]] = 2.0 # Minor allele
    if alleles > 2:
        trans[count[2][0]] = 1.0 # Heterozygous
    trans[None] = None
    return [ trans[g] for g in gs ]
