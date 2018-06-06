#!/usr/bin/env python

#
# 2015-09-02
# Ewa A. Bergmann (ewa.a.bergmann@gmail.com)
# New York Genome Center


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
import numpy as np
from collections import defaultdict
from MathOperations import *
from math import log10
import scipy
from scipy import stats


baseQ_max = 60


def create_conditional_likelihood_of_base_dict(checkpoints):
    
    D = defaultdict(lambda: defaultdict(dict))
    for bq in range(0, baseQ_max+1):
        Q = phred2prob(bq)/3

        f_AAAA_A = lambda x: (1-Q)
        f_AAAA_B = lambda x: Q
        f_BBBB_B = lambda x: (1-Q)
        f_ABAB_A = lambda x: 0.5
        f_ABAB_B = lambda x: 0.5
        
        f_AABB_A = lambda x: ((1-x) * (1 - Q)) + (x * Q)
        f_AABB_B = lambda x: (x * (1 - Q)) + ((1-x) * Q)
        f_AABA_A = lambda x: (1-0.5*x) * (1 - Q) + (0.5*x * Q)
        f_AABA_B = lambda x: (0.5*x) * (1 - Q) + (1-0.5*x) * Q
            
        f_ABAA_A = lambda x: (0.5 + 0.5*x) * (1 - Q) + (1-x)*0.5*Q
        f_ABAA_B = lambda x: (1-x)*0.5 * (1 - Q) + (0.5 + 0.5*x)  * Q
        for v in checkpoints:
            D['AABB_A'][v][bq] = log10(np.float64(f_AABB_A(v)))
            D['AABB_B'][v][bq] = log10(np.float64(f_AABB_B(v)))
            
            D['AABA_A'][v][bq] = log10(np.float64(f_AABA_A(v)))
            D['AABA_B'][v][bq] = log10(np.float64(f_AABA_B(v)))
            
            D['ABAA_A'][v][bq] = log10(np.float64(f_ABAA_A(v)))
            D['ABAA_B'][v][bq] = log10(np.float64(f_ABAA_B(v)))
            
            D['AAAA_A'][v][bq] = log10(np.float64(f_AAAA_A(v)))
            D['AAAA_B'][v][bq] = log10(np.float64(f_AAAA_B(v)))
            
            D['ABAB_A'][v][bq] = log10(np.float64(f_ABAB_A(v)))
            D['ABAB_B'][v][bq] = log10(np.float64(f_ABAB_B(v)))
            
    return(D)



def likelihood_per_marker(A, Scores, checkpoints, ref_basequals, alt_basequals):
    
    for bq in ref_basequals:
        for i in range(0, len(checkpoints)):
            v = checkpoints[i]
            ref_part = [Scores['AAAA_A'][v][bq], Scores['AABB_A'][v][bq], Scores['AABA_A'][v][bq], Scores['ABAB_A'][v][bq], Scores['ABAA_A'][v][bq], Scores['ABAA_B'][v][bq], Scores['AAAA_B'][v][bq], Scores['AABB_B'][v][bq], Scores['AABA_B'][v][bq]]
            A[i] += ref_part
    for bq in alt_basequals:
        for i in range(0, len(checkpoints)):
            v = checkpoints[i]
            alt_part = [Scores['AAAA_B'][v][bq], Scores['AABB_B'][v][bq], Scores['AABA_B'][v][bq], Scores['ABAB_B'][v][bq], Scores['ABAA_B'][v][bq], Scores['ABAA_A'][v][bq], Scores['AAAA_A'][v][bq], Scores['AABB_A'][v][bq], Scores['AABA_A'][v][bq]]
            A[i] += alt_part
    return(A)


def calculate_contamination_likelihood(checkpoints, Data, Scores):
    
    D = np.zeros([len(checkpoints), 1],dtype='float64')
    for marker_data in Data:
        A = np.zeros([len(checkpoints), 9],dtype='float64')
        for i in range(0, len(checkpoints)):
            A[i] += marker_data[0]
        A = likelihood_per_marker(A, Scores, checkpoints, marker_data[1], marker_data[2])
        for i in range(0, len(checkpoints)):
            D[i] += log10p(sum([pow(10,x ) for x in A[i]]))
            
    return(D)

def apply_brents_algorithm(Data, Scores, x1, x2, x3):
    def f(x):
        if type(x) is np.ndarray:
            x = x[0]
        ##print x
        Scores = create_conditional_likelihood_of_base_dict([x])
        D = calculate_contamination_likelihood([x], Data, Scores)
        return(D[0]*-1)
    
    if x1 == 0.0 and f(x2) > f(x1):
        return(x1)
    elif x3 == 1.0 and f(x2) > f(x3):
        return(x3)
        
    optimal_val = scipy.optimize.brent(f, brack=(x1, x2, x3), tol=1.0e-07, maxiter=30)
    if type(optimal_val) is np.ndarray:
        optimal_val = optimal_val[0]
    return(optimal_val)
    

