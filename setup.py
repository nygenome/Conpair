#!/usr/bin/env python

import glob
import os
from os.path import join, relpath
from setuptools import setup, find_packages

version = '0.2'


def find_package_files(dirpath, package, skip_exts=None):
    paths = []
    for (path, dirs, fnames) in os.walk(join(package, dirpath)):
        for fname in fnames:
            if skip_exts and any(fname.endswith(ext) for ext in skip_exts):
                continue
            fpath = join(path, fname)
            paths.append(relpath(fpath, package))
    return paths


setup(
    name='conpair',
    version=version,
    author='Ewa A Grabowska',
    description="Conpair: concordance and contamination estimator for tumor-normal pairs",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    keywords='bioinformatics',
    url='https://github.com/nygenome/Conpair',
    download_url='https://github.com/nygenome/Conpair/releases',
    packages=['conpair'],
    include_package_data=True,
    # data_files=[
    #     ('conpair', glob.glob('data/markers/*[.bed,.txt]'))
    # ],
    package_data={
        'conpair': find_package_files('', 'conpair')
    },
    zip_safe=False,
    install_requires=[
        'numpy>=1.7.0 ',
        'scipy>=0.14.0 ',
    ],
    scripts=[
        'scripts/estimate_tumor_normal_contamination.py',
        'scripts/run_gatk_pileup_for_sample.py',
        'scripts/verify_concordance.py',
    ],
    classifiers=[
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
