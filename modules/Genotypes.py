#!/usr/bin/env python

#
# 2015-10-08
# Ewa A. Bergmann (ewa.a.bergmann@gmail.com)
# New York Genome Center
#

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A. Bergmann (ewa.a.bergmann@gmail.com)


import sys
import os
from collections import defaultdict
import numpy as np
import random


def RAF2genotypeProb(RAF):
    ALT = 1 - RAF
    p_refref = RAF*RAF
    p_refalt = 2*RAF*ALT
    p_altalt = ALT*ALT
    return(p_refref, p_refalt, p_altalt)


def downsample(baseqs):
    '''Randomly downsample so that floats do not get reduced to 0.0'''
    if len(baseqs) > 450:
        x = float(450) / len(baseqs)
        return random.sample(baseqs, int(x*len(baseqs)))
    return baseqs


def compute_genotype_likelihood(ref_baseq, alt_baseq, normalize=True):
    '''Randomly downsample to 450x if needed and then get likelihood from base quals'''
    ref_baseq = downsample(ref_baseq)
    alt_baseq = downsample(alt_baseq)
    AA = 1
    BB = 1
    for i in ref_baseq:
        p = phred_to_p(i)
        AA *= (1 - p)
        BB *= p
    for i in alt_baseq:
        p = phred_to_p(i)
        AA *= p
        BB *= (1 - p)
        
    AB = 0.5**(len(ref_baseq) + len(alt_baseq))

    if normalize is True:
        S = AA + BB + AB
        AA = AA/S
        AB = AB/S
        BB = BB/S
    return(AA, AB, BB)

def phred_to_p(phred):
    return(np.float64(10**(-float(phred)/10)))


def prior_genotype_probability(RAF):
    return(RAF**2, 2*RAF*(1-RAF), (1-RAF)**2)
