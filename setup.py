#!/usr/bin/env python

from setuptools import setup
from pip.req import parse_requirements
from pip.download import PipSession
import sys

description = ("This library provides functionality in parsing and assessing the "
		"equivalency of different variant names according to recommendations by "
		"the Human Genome Variation Society (HGVS).")

def main():

    setup(
        name='hgvs-lib',
        version='1.0.0',
        description='Parses and compares two HGVS variant descriptions',
    	license='Creative Commons Attribution-Noncommercial-Share Alike license',
        long_description=open('README.txt').read(),
        author='Jennifer Yen',
        author_email='jennifer.yen@personalis.com',
        packages=['hgvs', 'tests'],
        include_package_data=True,
        package_data={
            '': ['requirements.txt'],
        },
        scripts=[],
        install_requires=['pip>=1.2'],
        tests_require=[str(line.req) for line in
                       parse_requirements('requirements.txt',
                                          session=PipSession())],
	classifiers = [
            'Development Status :: 4 - Beta',
	    'Intended Audience :: Developers',
	    'Intended Audience :: Science/Research',
	    'Programming Language :: Python',
	    'Programming Language :: Python :: 2',
	    'Programming Language :: Python :: 2.6',
	    'Programming Language :: Python :: 2.7',
	    'Programming Language :: Python :: 3',
	    'Programming Language :: Python :: 3.2',
	    'Programming Language :: Python :: 3.3',
	    'Programming Language :: Python :: 3.4',
	    'Topic :: Scientific/Engineering :: Bio-Informatics',
	    ],
	keywords='bioinformatics'
     )

if __name__ == '__main__':
    main()
