#!/usr/bin/env python2.7

#
# Version 0.1 2015-09-02
# Ewa Grabowska (egrabowska@nygenome.org)
# New York Genome Center
#

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2014) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Grabowska (egrabowska@nygenome.org)


import os
import itertools
from Genotypes import *
from collections import OrderedDict


class Marker:
    
    def __init__(self, chrom, pos, ref, alt, RAF):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.RAF = RAF


def get_markers(marker_file):
    Marker_Set = OrderedDict()
    if not os.path.isfile(marker_file):
        sys.stderr.write("The specified marker file cannot be find: " + marker_file)
    markers = open(marker_file)
    for line in markers:
        if "#" in line:
            continue
        line = line.split()
        M = Marker(line[0], line[1], line[2], line[3], float(line[4]))
        Marker_Set[line[0] + ":" + line[1]] = M
    markers.close()
    return(Marker_Set)


class Pileup:
    
    def __init__(self, chrom, pos, ref, qual_A, qual_C, qual_G, qual_T):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.Quals = {}
        self.Quals['A'] = qual_A
        self.Quals['C'] = qual_C
        self.Quals['G'] = qual_G
        self.Quals['T'] = qual_T
        self.depth = sum([len(qual_A), len(qual_C), len(qual_G), len(qual_T)])


def parse_mpileup_line(line, min_map_quality=0, min_base_quality=0):
    
    line = line.split()
    chrom = line[0]
    pos = line[1]
    ref = line[2]
    bases = line[3]
    baseQs = baseQ2int(line[4])
    
    if min_map_quality > 0:
        verbose_lines = line[6].split(',')
        mapqs = [v.split('@')[-1] for v in verbose_lines]
        mapqs_above_threshold = set([i for i in xrange(0, len(mapqs)) if int(mapqs[i]) >= min_map_quality])

    indexes_A = find_all_positions_of_char(bases, 'A')
    indexes_C = find_all_positions_of_char(bases, 'C')
    indexes_G = find_all_positions_of_char(bases, 'G')
    indexes_T = find_all_positions_of_char(bases, 'T')
    
    if min_map_quality > 0:
        indexes_A = mapqs_above_threshold.intersection(indexes_A)
        indexes_C = mapqs_above_threshold.intersection(indexes_C)
        indexes_G = mapqs_above_threshold.intersection(indexes_G)
        indexes_T = mapqs_above_threshold.intersection(indexes_T)
    
    A_base_quals_list = [baseQs[i] for i in indexes_A if int(baseQs[i]) >= min_base_quality]
    C_base_quals_list = [baseQs[i] for i in indexes_C if int(baseQs[i]) >= min_base_quality]
    G_base_quals_list = [baseQs[i] for i in indexes_G if int(baseQs[i]) >= min_base_quality]
    T_base_quals_list = [baseQs[i] for i in indexes_T if int(baseQs[i]) >= min_base_quality]
    
    P = Pileup(chrom, pos, ref, A_base_quals_list, C_base_quals_list, G_base_quals_list, T_base_quals_list)
    return(P)


def genotype_likelihoods_for_markers(Markers, mpileup_file, min_map_quality=0, min_base_quality=0):
    
    M = dict()
    f = open(mpileup_file)
    
    for line in f:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = parse_mpileup_line(line, min_map_quality=min_map_quality, min_base_quality=min_base_quality)
        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue
    
        ref = marker.ref
        alt = marker.alt
        RAF = marker.RAF
        
        if pileup.Quals[ref] == [] and pileup.Quals[alt] == []:
            M[pileup.chrom + ":" + pileup.pos] = None
            continue
        
        AA_likelihood, AB_likelihood, BB_likelihood = compute_genotype_likelihood(pileup.Quals[ref], pileup.Quals[alt], normalize=False)
        prAA, prAB, prBB = prior_genotype_probability(RAF)
        
        
        M[pileup.chrom + ":" + pileup.pos] = {'likelihoods' : [AA_likelihood/prAA, AB_likelihood/prAB, BB_likelihood/prBB], 'coverage': pileup.depth}
        
    for m in Markers:
        try:
            v = M[m]
        except:
            M[m] = None
    
    f.close()
    return(M)
    
    



def pileup2acgt(pileup, ref):
    ref = ref.upper()
    nts = ''
    i = 0
    while i < len(pileup):
        r = pileup[i]
        if r in ['$',"!",']','>']:
            i += 1
            continue
        if r == "^":
            i += 2
            continue
        if r == "+" or r == "-":
            indel_len = int(''.join(itertools.takewhile(str.isdigit, pileup[i+1:])))
            i += 1 + len(str(indel_len)) + indel_len
            continue
        if r == ',' or r == '.':
            r = ref
        else:
            r = r.upper()
            if r not in ['A', 'C', 'G', 'T', '*', 'N']:
                i += 1
                continue
        nts += r
        i += 1
    return(nts)


def baseQ2int(baseQ_string, scaling_factor=33):
    ints = []
    for bq in baseQ_string:
        bq = ord(bq) - scaling_factor
        ints.append(bq)
    return(ints)


def find_all_positions_of_char(s, char):
    indexes = [i for i in xrange(0, len(s)) if s[i] == char]
    return(indexes)





