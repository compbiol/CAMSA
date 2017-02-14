from setuptools import setup
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import camsa

setup(
    name='CAMSA',
    version=camsa.VERSION,
    author='Sergey Aganezov',
    author_email='aganezov@gwu.edu',
    description='CAMSA: a tool for Comparative Analysis and Merging of Scaffold Assemblies',
    license='MIT',
    keywords="comparative genomics, scaffolding, genome assembly",
    url='https://github.com/compbiol/camsa',
    packages=['', 'camsa', 'camsa.core', 'camsa.utils', 'camsa.utils.camsa', 'camsa.utils.fasta', 'camsa.utils.agp', 'camsa.utils.ragout', 'camsa.utils.grimm'],
    include_package_data=True,
    install_requires=['six>=1.10.0', 'networkx>=1.11', 'Jinja2>=2.8', 'enum34>=1.1.6', 'blist>=1.3.6', 'ConfigArgParse>=0.10.0',
                      'biopython>=1.67', 'bg>=1.8.1'],
    scripts=["camsa/run_camsa.py",
             "camsa/utils/ragout/ragout_coords2fasta.py", "camsa/utils/ragout/ragout_coords_coverage.py", "camsa/utils/ragout/ragout_coords2camsa_seqi.py", "camsa/utils/ragout/ragout_coords2camsa_points.py",
             "camsa/utils/grimm/grimm2camsa_points.py",
             "camsa/utils/fasta/fasta2camsa_points.py", "camsa/utils/fasta/fasta2camsa_seqi.py", "camsa/utils/fasta/camsa_points2fasta.py",
             "camsa/utils/agp/agp2camsa_points.py"],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    long_description=open("pypi_full.rst").read()

)
