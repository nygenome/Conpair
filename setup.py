#!/usr/bin/env python
"""
MultiQC_bcbio is a plugin for MultiQC, providing additional tools which are
specific to bcbio_nextgen.py pipeline.

For more information about bcbio_nextgen, see http://github.com/chapmanb/bcbio_nextgen
For more information about MultiQC, see http://multiqc.info
"""
import glob
from setuptools import setup, find_packages

version = '0.2'

setup(
    name = 'conpair',
    version = version,
    author = 'Ewa A Grabowska',
    description = "Conpair: concordance and contamination estimator for tumor-normal pairs",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    keywords = 'bioinformatics',
    url = 'https://github.com/nygenome/Conpair',
    download_url = 'https://github.com/nygenome/Conpair/releases',
    packages = ['conpair'],
    include_package_data = True,
    data_files = [
        ('conpair', glob.glob('data/markers/*[.bed,.txt]'))
    ],
    zip_safe=False,
    install_requires = [
        'numpy>=1.7.0 ',
        'scipy>=0.14.0 ',
    ],
    scripts = [
        'scripts/estimate_tumor_normal_contamination.py',
        'scripts/run_gatk_pileup_for_sample.py',
        'scripts/verify_concordance.py',
    ],
    classifiers = [
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)

