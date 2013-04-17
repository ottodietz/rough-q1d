# -*- coding: utf-8 -*-
from distutils.core import setup

setup(
    name='rough-q1d',
    version='0.1',
    packages=['q1d',],
    license='GPLv3'
    long_description=open('README.md').read(),
    author=['Jörg Doppler, Otto Dietz']
    author_email=['joerg.doppler@tuwien.ac.at','otto.dietz@physik.hu-berlin.de']
    packages=['q1d']
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='https://github.com/ottodietz/rough-q1d'
    license='LICENSE.txt'
    description='Tools to analyse rough boundary in quasi-one dimensional systems for python.'
    long_description=open('README.md').read(),
    install_requires=[ "numpy", "scipy" ]
)
