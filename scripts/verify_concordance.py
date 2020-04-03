#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.15
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

import sys
import os
import optparse
import math
from collections import defaultdict
from conpair.ContaminationModel import *
from conpair.ContaminationMarker import *
from conpair import find_markers_file


desc = """Program to verify tumor-normal sample concordance"""
parser = optparse.OptionParser(version='%prog version 0.15 3/August/2016', description=desc)
parser.add_option('-T', '--tumor_pileup', help='TUMOR PILEUP FILE [mandatory field]', action='store')
parser.add_option('-N', '--normal_pileup', help='NORMAL PILEUP FILE [mandatory field]', action='store')
parser.add_option('-M', '--markers', help='MARKER FILE [Conpair-GRCh37-default]', action='store')
parser.add_option('-g', '--genome', help='Instead of the marker file path, you can specify the genome build name (GRCh37, GRCh38, GRCm38) [default: GRCh37]', default='GRCh37', action='store')
parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
parser.add_option('-O', '--outfile', help='TXT OUTPUT FILE [stdout by default]', default="-", type='string', action='store')
parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [default: 10]', default=10, type='int', action='store')
parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY [default: 20]', default=20, type='int', action='store')
parser.add_option('-H', '--normal_homozygous_markers_only', help='USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)', default=False, action='store_true')

(opts, args) = parser.parse_args()

if not opts.tumor_pileup or not opts.normal_pileup:
    parser.print_help()
    sys.exit(1)
    
if not os.path.exists(opts.tumor_pileup):
    print('ERROR: Input tumor file {0} cannot be found.'.format(opts.tumor_pileup))
    sys.exit(1)
    
if not os.path.exists(opts.normal_pileup):
    print('ERROR: Input normal file {0} cannot be found.'.format(opts.normal_pileup))
    sys.exit(1)

markers_file = find_markers_file(opts, '.txt')

Markers = get_markers(markers_file)
COVERAGE_THRESHOLD = opts.min_cov
MMQ = opts.min_mapping_quality
MBQ = opts.min_base_quality
AA_BB_only = opts.normal_homozygous_markers_only

Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, opts.normal_pileup, min_map_quality=MMQ, min_base_quality=MBQ)
Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, opts.tumor_pileup, min_map_quality=MMQ, min_base_quality=MBQ)

concordant = 0
discordant = 0
for m in Markers:
    NL = Normal_genotype_likelihoods[m]
    TL = Tumor_genotype_likelihoods[m]
    if NL is None or TL is None:
        continue
    if NL['coverage'] < COVERAGE_THRESHOLD or TL['coverage'] < COVERAGE_THRESHOLD:
        continue
    if AA_BB_only:
        if NL['likelihoods'].index(max(NL['likelihoods'])) == 1:
            continue
    if NL['likelihoods'].index(max(NL['likelihoods'])) == TL['likelihoods'].index(max(TL['likelihoods'])):
        concordant += 1
    else:
        discordant += 1

if concordant+discordant == 0:
    print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(COVERAGE_THRESHOLD))
    sys.exit(0)


if opts.outfile == "-":
    print(round(float(concordant)/(concordant+discordant), 3))
    print("Based on " + str(concordant+discordant) + "/" + str(len(Markers)) + " markers (coverage per marker threshold: " + str(COVERAGE_THRESHOLD) + " reads)")
    print("Minimum mappinq quality: " + str(MMQ))
    print("Minimum base quality: " + str(MBQ))
else:
    outfile = open(opts.outfile, 'w')
    outfile.write("Concordance: "+ str(round(100.0*float(concordant)/(concordant+discordant), 2)) + "%\n")
    outfile.write("Based on " + str(concordant+discordant) + "/" + str(len(Markers)) + " markers (coverage per marker threshold : " + str(COVERAGE_THRESHOLD) + " reads)\n")
    outfile.write("Minimum mappinq quality: " + str(MMQ) + "\n")
    outfile.write("Minimum base quality: " + str(MBQ) + "\n")
    outfile.close()
                  

