from setuptools import setup

setup(
    name='CAMSA',
    version='0.9.0',
    author='Sergey Aganezov',
    author_email='aganezov@gwu.edu',
    description='CAMSA: a tool for Comparative Analysis and Merging of Scaffold Assemblies',
    license='MIT',
    keywords="comparative genomics, scaffolding, genome assembly",
    url='https://github.com/aganezov/camsa',
    packages=['core', 'utils', 'utils.camsa', 'utils.fasta'],
    package_dir={'': 'camsa'},
    install_requires=['six>=1.10.0', 'networkx>=1.11', 'Jinja2>=2.8', 'enum34>=1.1.6', 'blist>=1.3.6', 'ConfigArgParse>=0.10.0'],
    extras_require={
        'fasta': ['biopython>=1.67'],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
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

)
