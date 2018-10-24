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
from conpair import which, find_markers_file


desc = """Program to run GATK Pileup on a single sample"""
parser = optparse.OptionParser(version='%prog version 1.0 21/March/2016', description=desc)
parser.add_option('-B', '--bam', help='BAMFILE [mandatory field]', action='store')
parser.add_option('-O', '--outfile', help='OUTPUT FILE (PILEUP) [mandatory field]', type='string', action='store')
parser.add_option('-D', '--conpair_dir', help='CONPAIR DIR [$CONPAIR_DIR by default]', action='store')
parser.add_option('-R', '--reference', help='REFERENCE GENOME [GRCh37 by default]', action='store')
parser.add_option('-M', '--markers', help='MARKER FILE [GRCh37-default]', action='store')
parser.add_option('-g', '--genome', help='Instead of the marker file path, you can specify the genome build name (GRCh37, GRCh38, GRCm38) [default: GRCh37]', default='GRCh37', action='store')
parser.add_option('-G', '--gatk', help='GATK JAR [by default, will check if gatk executable as available in PATH assuming it\'s GATK4 installed with bioconda. '
                                       'Otherwise, will assume GATK JAR is in $GATK_JAR environment variable]', action='store')
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

GATK = None
GATK4 = None
if opts.gatk:
    GATK = opts.gatk
elif 'GATK_JAR' in os.environ:
    GATK = os.environ['GATK_JAR']
elif which('gatk'):
    GATK4 = 'gatk'
else:
    print('ERROR: cannot find GATK. Try either installing GATK4 with conda (conda install -c bioconda gatk4), or provide GATK2 or GATK3 JAR with `--gatk` or GATK_JAR environment variable')
    sys.exit(2)

if GATK and not os.path.exists(GATK):
    print('ERROR: GATK jar {0} cannot be find.'.format(GATK))
    sys.exit(2)

conpair_dir = os.environ.get('CONPAIR_DIR', opts.conpair_dir)
markers_file = find_markers_file(opts, '.bed', conpair_dir=conpair_dir)

if opts.reference:
    reference_fa = opts.reference
    if not os.path.exists(reference_fa):
        print('ERROR: Reference genome {0} cannot be find.'.format(reference_fa))
        sys.exit(3)
else:
    if conpair_dir:
        reference_fa = os.path.join(conpair_dir, 'data', 'genomes', 'human_g1k_v37.fa')
        if not os.path.exists(reference_fa):
            print('ERROR: Please provide reference fasta file with --reference, or put it into ' + reference_fa)
            sys.exit(3)
    else:
        print('ERROR: Please provide reference fasta file with --reference')
        sys.exit(3)

if opts.temp_dir_java:
    JAVA_TEMP = "-Djava.io.tmpdir={}".format(opts.temp_dir_java)
    if not os.path.isdir(opts.temp_dir_java):
        os.makedirs(opts.temp_dir_java, 770)
else:
    JAVA_TEMP = ""

if GATK4:
    command_line = ("{GATK4} --java-options '{JAVA_TEMP} -Xmx{opts.xmx_java}' Pileup -R {reference_fa} -I {opts.bam} -L {markers_file} -O {opts.outfile} " +
                    "-verbose -RF NotDuplicateReadFilter -RF CigarContainsNoNOperator " +
                    "-RF MatchingBasesAndQualsReadFilter").format(**locals())
else:
    command_line = ("{opts.java} {JAVA_TEMP} -Xmx{opts.xmx_java} -jar {GATK} -T Pileup -R {reference_fa} -I {opts.bam} -L {markers_file} -o {opts.outfile} " +
                    "-verbose -rf DuplicateRead --filter_reads_with_N_cigar " +
                    "--filter_mismatching_base_and_quals").format(**locals())
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
