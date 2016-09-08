#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)


import os
import sys
import optparse
from shutil import move
from tempfile import NamedTemporaryFile


desc = """Program to run GATK Pileup on a single sample"""
parser = optparse.OptionParser(version='%prog version 1.0 21/March/2016', description=desc)
parser.add_option('-B', '--bam', help='BAMFILE [mandatory field]', action='store')
parser.add_option('-O', '--outfile', help='OUTPUT FILE (PILEUP) [mandatory field]', type='string', action='store')
parser.add_option('-D', '--conpair_dir', help='CONPAIR DIR [$CONPAIR_DIR by default]', action='store')
parser.add_option('-R', '--reference', help='REFERENCE GENOME [GRCh37 by default]', action='store')
parser.add_option('-M', '--markers', help='MARKER FILE [GRCh37-default]', action='store')
parser.add_option('-G', '--gatk', help='GATK JAR [$GATK by default]', action='store')
parser.add_option('-J', '--java', help='PATH to JAVA [java by default]', default='java', action='store')
parser.add_option('-t', '--temp_dir_java', help='temporary directory to set -Djava.io.tmpdir', action='store')
parser.add_option('-m', '--xmx_java', help='Xmx java memory setting [default: 12g]', default='12g', action='store')
parser.add_option('--remove_chr_prefix', help='REMOVE CHR PREFIX FROM THE CHROMOSOME COLUMN IN THE OUTPUT FILE [false by default]', default=False, action='store_true')


(opts, args) = parser.parse_args()

if not opts.bam or not opts.outfile:
    parser.print_help()
    sys.exit(1)

if not os.path.exists(opts.bam):
    print('ERROR: Specified bamfile {0} cannot be found.'.format(opts.bam))
    sys.exit(1)

if opts.conpair_dir:
    CONPAIR_DIR = opts.conpair_dir
else:
    CONPAIR_DIR = os.environ['CONPAIR_DIR']

if opts.gatk:
    GATK = opts.gatk
else:
    GATK = os.environ['GATK_JAR']

if not os.path.exists(GATK):
    print('ERROR: GATK jar {0} cannot be find.'.format(GATK))
    sys.exit(2)

if opts.markers:
    MARKER_FILE = opts.markers
else:
    MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed')

if not os.path.exists(MARKER_FILE):
    print('ERROR: Marker file {0} cannot be find.'.format(MARKER_FILE))
    sys.exit(2)

if opts.reference:
    REFERENCE = opts.reference
else:
    REFERENCE = os.path.join(CONPAIR_DIR, 'data', 'genomes', 'human_g1k_v37.fa')

if not os.path.exists(REFERENCE):
    print('ERROR: Reference genome {0} cannot be find.'.format(REFERENCE))
    sys.exit(3)

if opts.temp_dir_java:
    JAVA_TEMP = "-Djava.io.tmpdir={}".format(opts.temp_dir_java)
    if not os.path.isdir(JAVA_TEMP):
        os.makedirs(JAVA_TEMP, 0770)
else:
    JAVA_TEMP = ""

command_line = "java {0} -Xmx{1} -jar {2} -T Pileup -R {3} -I {4} -L {5} -o {6} -verbose -rf DuplicateRead --filter_reads_with_N_cigar --filter_mismatching_base_and_quals".format(JAVA_TEMP, opts.xmx_java, GATK, REFERENCE, opts.bam, MARKER_FILE, opts.outfile)
os.system(command_line)


if opts.remove_chr_prefix:
    print("Removing 'chr' prefix...")

    with NamedTemporaryFile(delete=False) as tmp_source:
        with open(opts.outfile) as source_file:
            for line in source_file:
                if line.startswith("chr"):
                    tmp_source.write(line[3:])
            	else:
                    tmp_source.write(line)

    move(tmp_source.name, source_file.name)
