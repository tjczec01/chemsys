# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 03:02:05 2020

@author: tjcze
"""

from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='symbolgen',
    url='https://github.com/tjczec01/symbolgen',
    author='Travis Czechorski',
    author_email='tjczec01@gmail.com',
    # Needed to actually package something
    packages=['symbolgen'],
    # Needed for dependencies
    install_requires=['numpy', 'sympy', 'IPython', 'tkinter', 'scipy'],
    # *strongly* suggested for sharing
    version='1.0',
    # The license can be anything you like
    license='MIT',
    description='An example of a python package from pre-existing code',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)