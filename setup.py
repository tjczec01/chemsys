# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 18:25:39 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

"""
import os
from setuptools import setup, find_packages
# from sphinx.setup_command import BuildDoc

# cmdclass = {'build_sphinx': BuildDoc}
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))

description_chemsys = str("""Interactive GUI based program that generates the overall species balance,
                  system of ODEs needed for the solve_ivp and odeint method,
                  and calculates the Jacobian both symbolically and numerically.
                  The resulting code can easily be copied and pasted as is to be integrated with the aforementioned SciPy functions.""")


with open(r"{}\README.md".format(cwd), "r") as fh:
    long_description = fh.read()

# with open(r'{}\LICENSE.txt'.format(cwd)) as f:
#     license = f.read()


setup(
    name="chemsys", # Replace with your own username
    version="1.0.43",
    author="Travis Czechorski",
    author_email="tjczec01@gmail.com",
    description="{}".format(description_chemsys),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=r"https://github.com/tjczec01/chemsys",
    packages = find_packages(),
    keywords = ['chemical engineering', 'chemistry', 'engineering',
                'chemical reactions', 'jacobian', 'ode', "Plug", "Flow", "Reactor", "PFR"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    python_requires='>=3.6')
