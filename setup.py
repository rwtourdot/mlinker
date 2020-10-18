#!/usr/bin/env python3
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='mlinker',
    description='mlinker',
    author='Rick Tourdot',
    url='https://github.com/rwtourdot/mlinker',
    download_url='https://github.com/rwtourdot/mlinker',
    author_email='rwtourdot@gmail.com',
    version='0.0.0',
    install_requires=[
        'pysam',
        'numpy',
        'pandas',
        'matplotlib',
        'PyVCF'
    ],
    scripts=[
        'scripts/compare_hic_10x_haplotype_single_chr.py',
        'scripts/energy_statistics.py',
        'scripts/extract_coverage_from_vcf.py',
        'scripts/remove_centromere_variants_vcf.py',
        'scripts/remove_low_quality_spin_variants.py'
    ],
)
