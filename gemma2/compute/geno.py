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
