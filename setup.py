#!/usr/bin/env python

import glob
import os

from setuptools import setup, find_packages

version = '0.1'

required = open('requirements.txt').read().split('\n')
with open(os.path.join(os.path.dirname(__file__), 'requirements.txt')) as f:
    install_requires = [line.rstrip() for line in f]


setup(
    name='seqbot',
    version=version,
    description='Scripts for sequencing automation',
    author='James Webber',
    author_email='james.webber@czbiohub.org',
    url='https://github.com/czbiohub/seqbot',
    packages=find_packages(),
    install_requires=install_requires,
    long_description='See https://github.com/czbiohub/seqbot',
    license=open("LICENSE").readline().strip(),
)
