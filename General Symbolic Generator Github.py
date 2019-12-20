# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 03:43:11 2019

@author: tjcze
"""

import sympy as sp
import numpy as np
import mpmath as mp
from mpmath import mpf
import matplotlib.pyplot as plt 
from sympy import plot_implicit, init_printing, preview, symbols, exp, init_session, simplify, latex, Function, derive_by_array, Add, Float, Mul, Integer, Pow, Symbol
from sympy.parsing.sympy_parser import parse_expr
from operator import mul
from functools import reduce
import decimal 
from decimal import Decimal, getcontext
import IPython.display as ip
from IPython.display import display, DisplayObject, display_latex, Math, Image, Latex
import pprint as pp
from PIL import Image, ImageDraw, ImageMath, ImageFont
from IPython.lib.latextools import latex_to_png
import matplotlib.patches as patches
from matplotlib import rcParams
import pyglet

sp.init_session(use_latex=True,quiet=True)

# chemicals = int(input("Enter number of chemical species --> "))
# chemical_names = []
# for i in range(0,chemicals,1):
#        name = str(input("Enter chemical {} name --> ".format(i+1)))
#        chemical_names.append(name)
# print(chemical_names)       
# reactions = int(input("Enter number of chemical reactions --> "))


Initreactions = [
{"Ea" : 1, "K_Value": 1,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 2, "K_Value": 2,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 3, "K_Value": 3,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 4, "K_Value": 4,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 1, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 5, "K_Value": 5,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}},
{"Ea" : 6, "K_Value": 6,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 1, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 7, "K_Value": 7,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}},
{"Ea" : 8, "K_Value": 8,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 1, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 9, "K_Value": 9,
 "Reactants" : {'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 10, "K_Value": 10,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 1, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 11, "K_Value": 11,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 1, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'VCM':  0}},
{"Ea" : 12, "K_Value": 12,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 13, "K_Value": 13,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  0}},
{"Ea" : 14, "K_Value": 14,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  0}},
{"Ea" : 15, "K_Value": 15,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  1, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 16, "K_Value": 16,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  1,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 17, "K_Value": 17,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  1}},
{"Ea" : 18, "K_Value": 18,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  1,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 19, "K_Value": 19,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 1, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 20, "K_Value": 20,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  2,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  1, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
{"Ea" : 21, "K_Value": 21,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  1,'C11' : 0, 'C112' : 0, 'R1': 2, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  2, 'Coke': 2,'CP':  0,'Di': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0,'C11' : 0, 'C112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'VCM':  0}},
]
Eqlist = [
{"Name" : 'EDC' ,
 "Reactions" : {"Reaction 1" : -1 ,"Reaction 2" : -1,"Reaction 3" : -1,"Reaction 4" : -1,"Reaction 5" : -1,"Reaction 6" : -1, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}}, 
{"Name": 'EC',
  "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 1,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : -1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 1, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
  "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'HCl',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 1,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : 1, "Reaction 8" : 1, "Reaction 9" : 1, "Reaction 10" : 1, "Reaction 11" :1, "Reaction 12" :0, "Reaction 13" :1, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" : 1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'Coke',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" : 1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'CP',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" : 1, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'Di',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" : 1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" : 1, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :-1, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C4H6Cl2',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" : 1, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C6H6',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" : 1/4, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C2H2',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :1, "Reaction 19" :0, "Reaction 20" : -1, "Reaction 21" :-1/2 }, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :-1, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C11',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : -1, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C112',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 0, "Reaction 11" :-1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R1',
 "Reactions" : {"Reaction 1" : 1 ,"Reaction 2" : -1,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : -1, "Reaction 8" : -1, "Reaction 9" : -1, "Reaction 10" : -1, "Reaction 11" : -1, "Reaction 12" : -1, "Reaction 13" : -1, "Reaction 14" : 0, "Reaction 15" : 1, "Reaction 16" : 1, "Reaction 17" :1, "Reaction 18" :1, "Reaction 19" : 1, "Reaction 20" : 1/4, "Reaction 21" : -1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :-1, "Reaction 18" :-1, "Reaction 19" :-1, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R2',
 "Reactions" : {"Reaction 1" : 1 ,"Reaction 2" : 0,"Reaction 3" : -1,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : -1, "Reaction 8" :0, "Reaction 9" : 1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : -1, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R3',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 1,"Reaction 3" : 1,"Reaction 4" : 1,"Reaction 5" : 1,"Reaction 6" : 1, "Reaction 7" :0, "Reaction 8" : -1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" : -1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R4',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : -1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 1, "Reaction 11" :0, "Reaction 12" : 1, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" : -1, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R5',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : -1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" : 1, "Reaction 14" : 1, "Reaction 15" :0, "Reaction 16" :-1, "Reaction 17" :0, "Reaction 18" :-1, "Reaction 19" :0, "Reaction 20" :-1/4, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :1, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R6',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : -1, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" : 1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :-1, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :1, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'VCM',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :1, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" : -1, "Reaction 13" : -1, "Reaction 14" : -1, "Reaction 15" : -1, "Reaction 16" : -1, "Reaction 17" : 1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :-1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}}, 
]

Initreactionsf = [
{"Ea" : 1, "K_Value": 1,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 2, "K_Value": 2,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  1, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 3, "K_Value": 3,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 4, "K_Value": 4,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 5, "K_Value": 5,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 1, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 6, "K_Value": 6,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}},
{"Ea" : 7, "K_Value": 7,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 1, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 8, "K_Value": 8,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  1,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 1, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 9, "K_Value": 9,
 "Reactants" : {'EDC':  1, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  1, 'VCM':  0}},
{"Ea" : 10, "K_Value": 10,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}},
{"Ea" : 11, "K_Value": 11,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 1, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 12, "K_Value": 12,
 "Reactants" : {'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 13, "K_Value": 13,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 1, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 14, "K_Value": 14,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 1, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 15, "K_Value": 15,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 1, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  1,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 16, "K_Value": 16,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  1, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 17, "K_Value": 17,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 18, "K_Value": 18,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  1, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 19, "K_Value": 19,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  1, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  1, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 20, "K_Value": 20,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  1, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 21, "K_Value": 21,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  1, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 22, "K_Value": 22,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  1, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  1}},
{"Ea" : 23, "K_Value": 23,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  1, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 24, "K_Value": 24,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 1, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 25, "K_Value": 25,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  1,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 1, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 26, "K_Value": 26,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  1, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  1, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 1, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 27, "K_Value": 27,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  1, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 1, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 28, "K_Value": 28,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'R7':  0,  'R8':  0, 'CCl4':  1, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 1, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 29, "K_Value": 29,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  1, 'R7':  0,  'R8':  1, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 1, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  1, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 30, "K_Value": 30,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  2, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  1, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  1, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 1, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
{"Ea" : 31, "K_Value": 31,
 "Reactants" : {'EDC':  0, 'EC':  0, 'HCl':  0, 'Coke': 0, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  1, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 2, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}, 
 "Products" :{'EDC':  0, 'EC':  0, 'HCl':  2, 'Coke': 2, 'CP':  0, 'Di': 0, 'Tri': 0, 'C4H6Cl2':  0, 'C6H6':  0, 'C2H2':  0, 'C11' : 0, 'C112' : 0, 'C1112' : 0, 'R1': 0, 'R2':  0, 'R3':  0, 'R4':  0, 'R5':  0, 'R6':  0, 'R7':  0,  'R8':  0, 'CCl4':  0, 'CHCl3':  0, 'VCM':  0}},
]
Eqlistf = [
{"Name" : 'EDC' , 
 "Reactions" : {"Reaction 1" : -1,"Reaction 2" : 0,"Reaction 3" : -1,"Reaction 4" : -1, "Reaction 5" : -1,"Reaction 6" : -1, "Reaction 7" : -1, "Reaction 8" : -1, "Reaction 9" : -1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'EC' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" : -1, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" : 1, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'HCl' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 1,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 1, "Reaction 11" : 1, "Reaction 12" : 1, "Reaction 13" : 1, "Reaction 14" : 1, "Reaction 15" : 1, "Reaction 16" : 1, "Reaction 17" :0, "Reaction 18" : 1, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" : 1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'Coke' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" : 1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'CP' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" : 1, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'Di' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" : 1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" : 1, "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" : 1, "Reaction 28" :0 , "Reaction 29" : 1, "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" : -1, "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'Tri' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" : 1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" : -1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C4H6Cl2' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" : 1, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C6H6' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" : 1, "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C2H2' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" : 1, "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" : -1, "Reaction 31" : -1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" : -1, "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C11' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" : -1, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C112' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : 1, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : -1, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" : 1, "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'C1112' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" : 1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" : -1, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : 1, "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : -1 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R1' ,
 "Reactions" : {"Reaction 1" : 1,"Reaction 2" : 1,"Reaction 3" : -1,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : -1, "Reaction 11" : -1, "Reaction 12" : -1, "Reaction 13" : -1, "Reaction 14" : -1, "Reaction 15" : -1, "Reaction 16" : -1, "Reaction 17" : -1, "Reaction 18" : -1, "Reaction 19" :0, "Reaction 20" : 1, "Reaction 21" : 1, "Reaction 22" : 1, "Reaction 23" : 1, "Reaction 24" : 1, "Reaction 25" : 1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" : 1, "Reaction 31" : -1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" : -1, "Reaction 23" : -1, "Reaction 24" : -1, "Reaction 25" : -1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R2' ,
 "Reactions" : {"Reaction 1" : 1,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : -1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : -1, "Reaction 11" :0, "Reaction 12" : 1, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" : -1, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R3' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 1,"Reaction 4" : 1,"Reaction 5" : 1,"Reaction 6" : 1, "Reaction 7" : 1, "Reaction 8" : 1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" : -1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" : -1, "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" : 1, "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R4' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : -1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" : 1, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" : 1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" : -1, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" : -1, "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R5' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" :  -1, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" : 1, "Reaction 19" : 1, "Reaction 20" :0, "Reaction 21" : -1, "Reaction 22" :0 , "Reaction 23" : -1, "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" : -1, "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" : -1, "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" : 1, "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R6' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : -1, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 1, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" : -1, "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : -1, "Reaction 29" : -1, "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" : 1, "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : 1, "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R7' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" : -1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" : 1, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" : -1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" : 1, "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'R8' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 1,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : -1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" : 1, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" : 1, "Reaction 27" : 1, "Reaction 28" : 1, "Reaction 29" : -1, "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : -1, "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'CCl4' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : -1,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" : -1, "Reaction 27" : -1, "Reaction 28" : -1, "Reaction 29" : 1, "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" : 1, "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'CHCl3' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" : -1, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" :0 , "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
{"Name" : 'VCM' ,
 "Reactions" : {"Reaction 1" : 0,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 1, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 1, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" : -1, "Reaction 18" : -1, "Reaction 19" : -1, "Reaction 20" : -1, "Reaction 21" : -1, "Reaction 22" : 1, "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0, "Reaction 22" : -1, "Reaction 23" :0 , "Reaction 24" :0 , "Reaction 25" :0 , "Reaction 26" :0 , "Reaction 27" :0 , "Reaction 28" :0 , "Reaction 29" :0 , "Reaction 30" :0 , "Reaction 31" :0}}, 
]

De = decimal.Decimal

namespd = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM','T0','T1','Pure']
namespdf = ['EDC','EC','HCl','Coke', 'CP','Di','Tri','C4H6Cl2','C6H6','C2H2','C11','C112','C1112','R1','R2','R3','R4','R5','R6','R7','R8','CCl4','CHCl3','VCM','T0','T1','Pure']

names = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM','T0','T1']
temps = ['0','1']
namese = ['EDC','EC','HCl','Coke', 'CP','Di','C_{4}H_{6}Cl_{2}','C_{6}H_{6}','C_{2}H_{2}','C11','C112','R_{1}','R_{2}','R_{3}','R_{4}','R_{5}','R_{6}','VCM','T_{0}','T_{1}']
namesef = ['EDC','EC','HCl','Coke', 'CP','Di','Tri','C_{4}H_{6}Cl_{2}','C_{6}H_{6}','C_{2}H_{2}','C11','C112','C1112','R_{1}','R_{2}','R_{3}','R_{4}','R_{5}','R_{6}','R_{7}','R_{8}','CCl_{4}','CHCl_{3}','VCM','T_{0}','T_{1}']
nameslatex = [sp.latex(sp.symbols(r'{}'.format(i))) for i in namese]
nameslatexf = [sp.latex(sp.symbols(r'{}'.format(i))) for i in namesef]
namesf = ['EDC','EC','HCl','Coke', 'CP','Di','Tri','C4H6Cl2','C6H6','C2H2','C11','C112','C1112','R1','R2','R3','R4','R5','R6','R7','R8','CCl4','CHCl3','VCM','T','dT/dz']
tempsf = ['0','1']
consf = ['1','2','3']

con1 = sp.symbols('con1')
con2 = sp.symbols('con2')
con3 = sp.symbols('con3')
T1 = sp.symbols('T1')
T0 = sp.symbols('T0')
consv = [5.228861948970217258E06, 4.573981356040483661E04, 2.392702360656382751E11]

Ea = [342.0,7.0,42.0,45.0,34.0,48.0,13.0,12.0,4.0,6.0,15.0,0.0,56.0,31.0,30.0,61.0,84.0,90.0,70.0,20.0,70.0] #[kJ/mol]
Eab = [De(x*1000) for x in Ea] #kJ/mol
k_0 = [5.9E15, 1.3E13, 1.0E12, 5.0E11, 1.2E13, 2.0E11, 1.0E13, 1.0E13, 1.7E13, 1.2E13, 1.7E13, 9.1E10, 1.2E14, 5.0E11, 2.0E10, 3.0E11, 2.1E14, 5.0E14, 2.0E13, 1.0E14, 1.6E14] #If first order: [1/s]. If second order: [cm^3/mol*s]
k_0b = [De(x) for x in k_0] 
k_0f = [5.90E15,	2.20E12,	1.30E13,	1.20E13,	1.00E12,	5.00E11,	2.00E11,	1.00E11,	1.00E12,	1.00E13,	1.00E13,	1.70E13,	1.20E13,	1.70E13,	1.70E13,	1.60E13,	9.10E10,	1.20E14,	3.00E11,	2.00E10,	5.00E11,	2.10E14,	5.00E14,	2.00E13,	2.50E13,	1.00E12, 5.00E11,	5.00E11,	1.00E13,	1.00E14,	1.60E14]
k_0fb = [De(x) for x in k_0f]
Eaf = [342,230,7,34,42,45,48,56,63,13,12,4,6,15,17,14,0,56,61,30,31,84,90,70,70,33,33,33,13,20,70]
Eafb = [De(x*1000) for x in Eaf]

Rvalc = 8.314 # [J/mol*K]
temp_c = 500 # [°C]
temp_k = temp_c + 273.15 # [°K]
Temp_K = De(temp_k)
Twallv = temp_c + 273.15
Twall = sp.symbols('Twall')
Rval = sp.symbols('R_val')
ksym = []
easym = []
elist = ['C_EDC',' C_EC',' C_HCl',' C_Coke',' C_CP',' C_Di',' C_Tri',' C_C4H6Cl2',' C_C6H6',' C_C2H2',' C_C11',' C_C112',' C_C1112',' C_R1',' C_R2',' C_R3',' C_R4',' C_R5',' C_R6',' C_R7',' C_R8',' C_CCl4',' C_CHCl3',' C_VCM',' T0',' T1']
elists=['C_EDC','C_EC','C_HCl','C_Coke','C_CP','C_Di','C_C4H6Cl2','C_C6H6','C_C2H2','C_C11','C_C112','C_R1','C_R2','C_R3','C_R4','C_R5','C_R6','C_VCM','T0','T1']
for ks in range(0,len(k_0),1):
    vv = sp.symbols('k_0s[{}]'.format(ks))
    ksym.append(vv)

for eg in range(0,len(Ea),1):
    ve = sp.symbols('Eas[{}]'.format(eg))
    easym.append(ve)

ksym[0], ksym[1], ksym[2], ksym[3], ksym[4], ksym[5], ksym[6], ksym[7], ksym[8], ksym[9], ksym[10], ksym[11], ksym[12], ksym[13], ksym[14], ksym[15], ksym[16], ksym[17], ksym[18], ksym[19], ksym[20] = sp.symbols('k_0[0] k_0[1] k_0[2] k_0[3] k_0[4] k_0[5] k_0[6] k_0[7] k_0[8] k_0[9] k_0[10] k_0[11] k_0[12] k_0[13] k_0[14] k_0[15] k_0[16] k_0[17] k_0[18] k_0[19] k_0[20]')
easym[0], easym[1], easym[2], easym[3], easym[4], easym[5], easym[6], easym[7], easym[8], easym[9], easym[10], easym[11], easym[12], easym[13], easym[14], easym[15], easym[16], easym[17], easym[18], easym[19], easym[20] = sp.symbols('Ea[0] Ea[1] Ea[2] Ea[3] Ea[4] Ea[5] Ea[6] Ea[7] Ea[8] Ea[9] Ea[10] Ea[11] Ea[12] Ea[13] Ea[14] Ea[15] Ea[16] Ea[17] Ea[18] Ea[19] Ea[20]')
C_EDC , C_EC , C_HCl , C_Coke , C_CP , C_Di , C_C4H6Cl2 , C_C6H6 , C_C2H2 , C_C11 , C_C112 , C_R1 , C_R2 , C_R3 , C_R4 , C_R5 , C_R6 , C_VCM , T0 , T1 = sp.symbols('C_EDC C_EC C_HCl C_Coke C_CP C_Di C_C4H6Cl2 C_C6H6 C_C2H2 C_C11 C_C112 C_R1 C_R2 C_R3 C_R4 C_R5 C_R6 C_VCM T0 T1')
eqn = sp.srepr(con1*T1 - (con2*(Twall - T0) +  con3*(ksym[0]*C_EDC*sp.exp((-1*sp.symbols('easym[0]'))/(Rval*T0)) + ksym[1]*C_EDC*C_R1*sp.exp(-easym[1]/(Rval*T0)) + ksym[2]*C_EDC*C_R2*sp.exp(-easym[2]/(Rval*T0)) + ksym[3]*C_EDC*C_R4*sp.exp(-easym[3]/(Rval*T0)) + ksym[4]*C_EDC*C_R5*sp.exp(-easym[4]/(Rval*T0)) + ksym[5]*C_EDC*C_R6*sp.exp(-easym[5]/(Rval*T0)))))
eqn2 = sp.srepr(con1*T1 - (con2*(Twall - T0) +  con3*(k_0[0]*C_EDC*sp.exp((-1*Ea[0])/(Rval*T0)) + k_0[1]*C_EDC*C_R1*sp.exp(-Ea[1]/(Rval*T0)) + k_0[2]*C_EDC*C_R2*sp.exp(-Ea[2]/(Rval*T0)) + k_0[3]*C_EDC*C_R4*sp.exp(-Ea[3]/(Rval*T0)) + k_0[4]*C_EDC*C_R5*sp.exp(-Ea[4]/(Rval*T0)) + k_0[5]*C_EDC*C_R6*sp.exp(-Ea[5]/(Rval*T0)))))

#print(eqn)
#print(eqn2)
def prod(seq):
    return reduce(mul, seq) if seq else 1

def symsgenerators(names,temps,cons):
    smmb = ["C_{}".format(i) for i in names]
    tsymbs = ["T_{}".format(j) for j in temps]
    csymbs = ["con{}".format(k) for k in cons]
    return smmb,tsymbs,csymbs

def symsgeneratorf(namesf,tempsf,cons):
    smmbf = ["C_{}".format(i) for i in namesf]
    tsymbsf = ["T_{}".format(j) for j in range(0,len(tempsf),1)]
    csymbsf = ["con{}".format(k) for k in cons]
    del smmbf[-1]
    del smmbf[-1]
#    smmbf.append(tsymbsf[0])
#    smmbf.append(tsymbsf[1])
#    lt = len(tsymbsf)
#    if lt !=0:
#        for qf,wf in enumerate(tsymbsf):
#            smmbf.append(tsymbsf[qf])
#        symsf = sp.symbols(smmbf, real=True, nonnegative=True)
#    else:
#        symsf = sp.symbols(smmbf, real=True, nonnegative=True)
    return smmbf,tsymbsf,csymbsf
def makesymbs(rxns,EQlist,names,temps,cons,consp):
    T1 = sp.symbols('T_1')
    T0 = sp.symbols('T_0')
    rvalsymb = sp.symbols('R_{val}')
    T0sym = sp.symbols('T_{0}')
    symbs = symsgenerators(names,temps,cons)[0]
    tsymbs = symsgenerators(names,temps,cons)[1]
    csymbs = symsgenerators(names,temps,cons)[2]
    EQs = []
    EQsf = []
    reacteqs = []
    prodeqs = []
    reacteqs2 = []
    prodeqs2 = []
    ks = []
    eas = []
    eqlistrc = []
    symbslist = list(symbs)
    for t in range(0,len(tsymbs),1):
        symbslist.append(tsymbs[t])
    for i,j in enumerate(Initreactions):
        productsb = []
        reactantsb = []
        Reactants = rxns[i]['Reactants']
        Products = rxns[i]['Products']
        eaval = rxns[i]['Ea'] - 1
        eaval2 = De(Eab[eaval-1])/(De(Rvalc)*T0)
        kk = rxns[i]['K_Value']-1
        kk2 = k_0f[kk]
        k_val = sp.symbols("k_0[{}]".format(rxns[i]['K_Value']-1))
        ea_val = sp.symbols("Ea_[{}]".format(eaval))
        ks.append(k_val)
        eas.append(ea_val)
        Eavaltop = ea_val
        Eavalbot = rvalsymb*T0sym
        rvals = list(Reactants.values())
        pvals = list(Products.values())
        for k,v  in zip(symbs,rvals): 
                v0 = Pow(Symbol(k),Integer(v))
                v1 = Mul(v0, Integer(v))
                if v1 != 0:
                    reactantsb.append(v1)
        for k2,v2  in zip(symbs,pvals):
                v0b = Pow(Symbol(k2),Integer(v2))
                v1b = Mul(v0b, Integer(v2))
                if v1b != 0:
                    productsb.append(v1b)
        llval = len(reactantsb)
        if llval == 2:
             finalr = sp.exp(Eavaltop/Eavalbot)*k_val*reactantsb[0] * reactantsb[1]
             finalr2 = sp.exp(eaval2)*De(kk2)*reactantsb[0] * reactantsb[1]
        elif llval == 1:
             finalr = sp.exp(Eavaltop/Eavalbot)*k_val*reactantsb[0]
             finalr2 = sp.exp(eaval2)*De(kk2)*reactantsb[0]
        reacteqs.append(finalr)
        reactantsb.clear()
        reacteqs2.append(finalr2)
        llval2 = len(productsb)
        if llval2 == 2:
             finalp =  sp.exp(Eavaltop/Eavalbot)*k_val*productsb[0]*productsb[1]
             finalp2 =  sp.exp(eaval2)*De(kk2)*productsb[0]*productsb[1]
        elif llval2 == 1:
             finalp =  sp.exp(Eavaltop/Eavalbot)*k_val*productsb[0]
             finalp2 =  sp.exp(eaval2)*De(kk2)*productsb[0]
        prodeqs.append(finalp)
        productsb.clear()
        prodeqs2.append(finalp2)
    for i,j in enumerate(Eqlist):
        eqlist = []
        eqlist2 = []
        rxn = Eqlist[i]["Reactions"]
        reverse = Eqlist[i]["Reverse"]
        nums3 = list(rxn.values())
        backs = list(reverse.values())
        for k,m in enumerate(nums3):
            forval = reacteqs[k]*m
            forval2 = reacteqs2[k]*m
            eqlist.append(forval)
            eqlist2.append(forval2)
        for h,g in enumerate(backs):
            reverseval = prodeqs[h]*g
            reverseval2 = prodeqs2[h]*g
            eqlist.append(reverseval)
            eqlist2.append(reverseval2)
        eqlistf = [sum(eqlist)]
        eqlistf2 = [sum(eqlist2)]
        eqlistrc.append(eqlist)
        eqlistr = [x*De(y) for x,y in zip(reacteqs2,nums3)]
        eqlistr2 = [x*De(y) for x,y in zip(prodeqs2,backs) if y != 0]
        lenb2 = len(eqlistr2)
        if lenb2 != 0:
            for i,j in enumerate(eqlistr2):
                eqlistr.append(eqlistr2[i])
        EQs.append(sp.sympify(eqlistf))
        EQsf.append(eqlistf2)
    def fmatrix(eqns):
        rval = len(eqns)
        cval = len(eqns[0])
        x = sp.MatrixSymbol('x', rval, cval)
        xlist = []
        X1 = sp.Matrix(x)
        M1 = sp.Matrix(eqns)
        for h in range(rval):
            for g in range(cval):
                vv1 = sp.symbols('x[{},{}]'.format(h,g))
                xlist.append(vv1)
        ff2 = lambda i,j: X1.xreplace({x:sp.sympify(M1)})
        for i in range(rval):
            for j in range(cval):
                    X2 = ff2(i,j)
        return X2
    M = fmatrix(eqlistrc)
    egslista = lambda l: [item for sublist in l for item in sublist]
    egslist = egslista(EQs)
    egslistb = lambda l2: [item2 for sublist2 in l2 for item2 in sublist2]
    egslist2 = egslistb(EQsf)
    var11 = symbslist
    prod2b = Mul(Symbol('T1'), Float(consp[0], precision=53))
    prod3b = Mul(Integer(-1), Float(consp[1], precision=53), Add(Mul(Integer(-1), Symbol('T_0')), Float(Twallv, precision=53)))
    prod4b = Mul(Integer(-1), Float(consp[2], precision=53), Add(Mul(Float(k_0[1], precision=53), Symbol('C_EDC'), Symbol('C_R1'), exp(Mul(Integer(-1), Float('0.84195333172961273', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Float(k_0[2], precision=53), Symbol('C_EDC'), Symbol('C_R2'), exp(Mul(Integer(-1), Float('5.0517199903776762', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Float(k_0[3], precision=53), Symbol('C_EDC'), Symbol('C_R4'), exp(Mul(Integer(-1), Float('5.4125571325475104', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Float(k_0[4], precision=53), Symbol('C_EDC'), Symbol('C_R5'), exp(Mul(Integer(-1), Float('4.0894876112581189', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Float(k_0[5], precision=53), Symbol('C_EDC'), Symbol('C_R6'), exp(Mul(Integer(-1), Float('5.7733942747173446', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Float(k_0[0], precision=53), Symbol('C_EDC'), exp(Mul(Integer(-1), Float('0.12027904738994467', precision=53), Float(Ea[0], precision=53), Pow(Symbol('T0'), Integer(-1)))))))
    prod5 = Add(prod3b,prod4b)
    prod6 = Add(prod2b, Mul(Integer(-1), prod5))
    eqfinal = prod6
    gg = Add(Mul(Symbol('T1'), Symbol('{}'.format(csymbs[0]))), Mul(Integer(-1), Symbol('{}'.format(csymbs[1])), Add(Mul(Integer(-1), Symbol('T0')), Symbol('Twall'))), Mul(Integer(-1), Symbol('{}'.format(csymbs[2])), Add(Mul(Symbol('C_EDC'), Symbol('C_R1'), Float(k_0[1]), exp(Mul(Integer(-1), Float(Ea[1]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Symbol('C_EDC'), Symbol('C_R2'), Float(k_0[2]), exp(Mul(Integer(-1), Float(Ea[2]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Symbol('C_EDC'), Symbol('C_R4'), Float(k_0[3]), exp(Mul(Integer(-1), Float(Ea[3]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Symbol('C_EDC'), Symbol('C_R5'), Float(k_0[4]), exp(Mul(Integer(-1), Float(Ea[4]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Symbol('C_EDC'), Symbol('C_R6'), Float(k_0[5]), exp(Mul(Integer(-1), Float(Ea[5]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Symbol('C_EDC'), Float(k_0[0]), exp(Mul(Integer(-1), Float(Ea[0]), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))))))
    EQs.append(gg)
    EQsf.append(eqfinal)
    egslist.append(T1)
    egslist2.append(T1)
    egslist.append(gg)
    egslist2.append(eqfinal)
    fj = sp.sympify(egslist)
    Jac = sp.zeros(len(egslist),len(var11))
    fj2 = sp.sympify(egslist2)
    Jac2 = sp.zeros(len(egslist2),len(var11))
    ksff = [*ks]
    for a, fi in enumerate(fj):
        for b, s in enumerate(var11):
            Jac[a,b] = sp.diff(fi, s)
    for a2, fi2 in enumerate(fj2):
        for b2, s2 in enumerate(var11):
            Jac2[a2,b2] = sp.diff(fi2, s2)
    

    return [egslist,egslist2, symbslist, ksff,eas,M,Jac,Jac2]


def makesymbsf(rxns,EQlist,namesf,tempsf,cons):
    Rval = sp.symbols('R_val')
    Twall = sp.symbols('T_wall')
    symbsf = symsgeneratorf(namesf,tempsf,cons)[0]
    symbsfnew = symbsf
    
    tsymsf = symsgeneratorf(namesf,tempsf,cons)[1]
    csymbsf = symsgeneratorf(namesf,tempsf,cons)[2]
    t0s = sp.symbols('{}'.format(tsymsf[0]))
    t1s = sp.symbols('{}'.format(tsymsf[1]))
    con1 = sp.symbols('{}'.format(csymbsf[0]))
    con2 = sp.symbols('{}'.format(csymbsf[1]))
    con3 = sp.symbols('{}'.format(csymbsf[2]))
    symbsfnew.append(t0s)
    symbsfnew.append(t1s)
    EQs = []
    EQsf = []
    reacteqs = []
    prodeqs = []
    reacteqs2 = []
    prodeqs2 = []
    ks = []
    eas = []
    eqlistrc = []
    for i,j in enumerate(Initreactionsf):
        productsb = []
        reactantsb = []
        Reactantsf = rxns[i]['Reactants']
        Productsf = rxns[i]['Products']
        eaval = rxns[i]['Ea'] - 1
        eaval2 = De(Eaf[eaval-1])/(De(Rvalc)*t0s)
        kk = rxns[i]['K_Value']-1
        kk2 = k_0f[kk]
        k_val = sp.symbols("k_0[{}]".format(rxns[i]['K_Value']-1))
        ea_val = sp.symbols("Ea_[{}]".format(eaval))
        ks.append(k_val)
        eas.append(ea_val)
        Eavaltop = ea_val
        Eavalbot = sp.combsimp(Rval*t0s)
        rvalsf = list(Reactantsf.values())
        pvalsf = list(Productsf.values())
        for k,v  in zip(symbsfnew,rvalsf): 
                v0 = Pow(Symbol(k),Integer(v))
                v1 = Mul(v0, Integer(v))
                if v1 != 0:
                    reactantsb.append(v1)
        for k2,v2  in zip(symbsfnew,pvalsf):
                v0b = Pow(Symbol(k2),Integer(v2))
                v1b = Mul(v0b, Integer(v2))
                if v1b != 0:
                    productsb.append(v1b)
        llval = len(reactantsb)
        if llval == 2:
             finalr = sp.exp(Eavaltop/Eavalbot)*k_val*reactantsb[0] * reactantsb[1]
             finalr2 = sp.exp(eaval2)*De(kk2)*reactantsb[0] * reactantsb[1]
        elif llval == 1:
             finalr = sp.exp(Eavaltop/Eavalbot)*k_val*reactantsb[0]
             finalr2 = sp.exp(eaval2)*De(kk2)*reactantsb[0]
        reacteqs.append(finalr)
        reactantsb.clear()
        reacteqs2.append(finalr2)
        llval2 = len(productsb)
        if llval2 == 2:
             finalp =  sp.exp(Eavaltop/Eavalbot)*k_val*productsb[0]*productsb[1]
             finalp2 =  sp.exp(eaval2)*De(kk2)*productsb[0]*productsb[1]
        elif llval2 == 1:
             finalp =  sp.exp(Eavaltop/Eavalbot)*k_val*productsb[0]
             finalp2 =  sp.exp(eaval2)*De(kk2)*productsb[0]
        prodeqs.append(finalp)
        productsb.clear()
        prodeqs2.append(finalp2)
    for i,j in enumerate(Eqlistf):
        eqlist = []
        eqlist2 = []
        rxn = Eqlistf[i]["Reactions"]
        reverse = Eqlistf[i]["Reverse"]
        nums3 = list(rxn.values())
        backs = list(reverse.values())
        for k,m in enumerate(nums3):
            forval = reacteqs[k]*m
            forval2 = reacteqs2[k]*m
            eqlist.append(forval)
            eqlist2.append(forval2)
        for h,g in enumerate(backs):
            reverseval = prodeqs[h]*g
            reverseval2 = prodeqs2[h]*g
            eqlist.append(reverseval)
            eqlist2.append(reverseval2)
        eqlistf = [sum(eqlist)]
        eqlistf2 = [sum(eqlist2)]
        eqlistrc.append(eqlist)
        eqlistr = [x*De(y) for x,y in zip(reacteqs2,nums3)]
        eqlistr2 = [x*De(y) for x,y in zip(prodeqs2,backs) if y != 0]
        lenb2 = len(eqlistr2)
        if lenb2 != 0:
            for i,j in enumerate(eqlistr2):
                eqlistr.append(eqlistr2[i])
        EQs.append(sp.sympify(eqlistf))
        EQsf.append(eqlistf2)
    def fmatrix(eqns):
        rval = len(eqns)
        cval = len(eqns[0])
        x = sp.MatrixSymbol('x', rval, cval)
        xlist = []
        X1 = sp.Matrix(x)
        M1 = sp.Matrix(eqns)
        for h in range(rval):
            for g in range(cval):
                vv1 = sp.symbols('x[{},{}]'.format(h,g))
                xlist.append(vv1)
        ff2 = lambda i,j: X1.xreplace({x:sp.sympify(M1)})
        for i in range(rval):
            for j in range(cval):
                    X2 = ff2(i,j)
        return X2
    M = fmatrix(eqlistrc)
    egslista = lambda l: [item for sublist in l for item in sublist]
    egslist = egslista(EQs)
    egslistb = lambda l2: [item2 for sublist2 in l2 for item2 in sublist2]
    egslist2 = egslistb(EQsf)
    var11f = symbsf
    var11f.append(t0s)
    var11f.append(t1s)
#    var11.append(csymbsf)
    prod0 = Add(Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R1', real=True, nonnegative=True), Symbol('k_0[2]'), exp(Mul(Symbol('Ea[2]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R2', real=True, nonnegative=True), Symbol('k_0[3]'), exp(Mul(Symbol('Ea[3]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R4', real=True, nonnegative=True), Symbol('k_0[4]'), exp(Mul(Symbol('Ea[4]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R5', real=True, nonnegative=True), Symbol('k_0[5]'), exp(Mul(Symbol('Ea[5]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R6', real=True, nonnegative=True), Symbol('k_0[6]'), exp(Mul(Symbol('Ea[6]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R7', real=True, nonnegative=True), Symbol('k_0[7]'), exp(Mul(Symbol('Ea[7]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R8', real=True, nonnegative=True), Symbol('k_0[8]'), exp(Mul(Symbol('Ea[8]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Symbol('C_EDC', real=True, nonnegative=True), Symbol('k_0[0]'), exp(Mul(Symbol('Ea[0]'), Pow(Symbol('Rval'), Integer(-1)), Pow(Symbol('T0'), Integer(-1))))))
    prod0f = Add(Mul(Integer(-1), Float('13000000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R1', real=True, nonnegative=True), exp(Mul(Float('27.664180899687274', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('12000000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R2', real=True, nonnegative=True), exp(Mul(Float('0.84195333172961273', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('1000000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R4', real=True, nonnegative=True), exp(Mul(Float('4.0894876112581189', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('500000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R5', real=True, nonnegative=True), exp(Mul(Float('5.0517199903776762', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('200000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R6', real=True, nonnegative=True), exp(Mul(Float('5.4125571325475104', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('100000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R7', real=True, nonnegative=True), exp(Mul(Float('5.7733942747173446', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('1000000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), Symbol('C_R8', real=True, nonnegative=True), exp(Mul(Float('6.7356266538369018', precision=53), Pow(Symbol('T0'), Integer(-1))))), Mul(Integer(-1), Float('5900000000000000.0', precision=53), Symbol('C_EDC', real=True, nonnegative=True), exp(Mul(Float('8.4195333172961266', precision=53), Pow(Symbol('T0'), Integer(-1))))))
    prod1 = sp.combsimp(con1*t1s) 
    prod2 = Add(Twall, Mul(Integer(-1), t0s))
    prod3 = sp.combsimp(con2*prod2) 
    prod4 = sp.combsimp(con3*prod0) 
    prod5 = Add(prod3,prod4)
    prod6 = Add(prod1, Mul(Integer(-1), prod5))
    eqfinal = prod6
    prod4f = sp.combsimp(con3*prod0f) 
    prod5f = Add(prod3,prod4f)
    prod6f = Add(prod1, Mul(Integer(-1), prod5f))
    eqfinalf = prod6f
    EQs.append(eqfinal)
    EQsf.append(eqfinalf)
    egslist.append(t1s)
    egslist2.append(t1s)
    egslist.append(eqfinal)
    egslist2.append(eqfinalf)
    fj = sp.sympify(egslist)
    Jac = sp.zeros(len(egslist),len(var11f))
    fj2 = sp.sympify(egslist2)
    Jac2 = sp.zeros(len(egslist2),len(var11f))
    ksff = [*ks]
    for a, fi in enumerate(fj):
        for b, s in enumerate(var11f):
            Jac[a,b] = sp.diff(fi, s)
    for a2, fi2 in enumerate(fj2):
        for b2, s2 in enumerate(var11f):
            Jac2[a2,b2] = sp.diff(fi2, s2)
    

    return [egslist,egslist2, symbsfnew, ksff,eas,M,Jac,Jac2]
#display(Latex(r'$\frac{d^2 T}{dz^2}$'))
ydot, ydotf, ys, kd, EA, IM, Jac, Jacf = makesymbs(Initreactions,Eqlist,names,temps,consf,consv)
ydot2, ydot2f, ys2, kd2, EA2, IM2, Jac2, Jac2f = makesymbsf(Initreactionsf,Eqlistf,namesf,tempsf,consf)  
startl = []
for i,j in enumerate(ys2):
       startl.append("{} = Decimal(C[{}])".format(j,i))
       print("{} = Decimal(C[{}])\n".format(j,i))
#print(startl)
#pp.pprint(sp.latex(ydot))
rcParams["text.usetex"] = True
latexeqns = []
for i,j in  enumerate(ydot):
       if i <= 17:
              indexk = i
              name = names[indexk]
              latexval = sp.latex(j)
              difname = str('dC_{}'.format(name))
              difz = str('dz')
              dxdz = str('{}/{}'.format(difname,difz))
              eqn = str('{}{}{}{}{}{}'.format('{',difname,'}','{','dz','}'))
              d4 = Latex(r'$\frac{}$'.format(eqn))
              d4b = Latex(r'${}$'.format(latexval))
       #       d5 = Latex(r'$\frac{}{}{}{}{}{}$ = '.format('{',difname,'}','{','dz','}') + r'${}$'.format(i))
              eqnval = str("{} = {}".format(sp.latex(sp.sympify(dxdz)),latexval))
              latexeqns.append(eqnval)
       elif i == 18:
              indexkf = i
              namef = namesf[i]
              latexval = sp.latex(j)
              difname = str('dT')
              difz = str('dz')
              dxdz = str('{}/{}'.format(difname,difz))
              eqn = str('{}{}{}{}{}{}'.format('{',difname,'}','{','dz','}'))
              d4 = Latex(r'$\frac{}$'.format(eqn))
              d4b = Latex(r'${}$'.format(latexval))
              eqnval = str("{} = {}".format(sp.latex(sp.sympify(dxdz)),latexval))
              latexeqns.append(eqnval)
       elif i == 19:
              indexkf = i
              namef = namesf[i]
              latexval = sp.latex(j)
              difname = str('d**2')
              difz = str('z**2')
              dxdz = str('{}T/d{}'.format(sp.latex(sp.sympify(difname)),sp.latex(sp.sympify(difz))))
              dxdz2 = str('{}T/d{}'.format(difname,difz))
#              print(dxdz2)
              eqn = str('{}{}{}{}{}{}'.format('{',difname,'}','{',difz,'}'))
              d4 = Latex(r'$\frac{}$'.format(eqn))
              d4b = Latex(r'${}$'.format(latexval))
              eqnval = str(r'\frac{d^2 T}{dz^2} = '+ r'{}'.format(latexval))
              latexeqns.append(eqnval)
#       print(eqn)
#       print(r'$\frac{}{}{}{}{}{}$ = '.format('{',difname,'}','{','dz','}') + r'${}$'.format(i))
latexeqnsf = []
for i2,j2 in enumerate(ydot2):
       if i2 <= 23:
              indexkf = i2
              namef = namesf[i2]
              latexvalf = sp.latex(j2)
              difnamef = str('dC_{}'.format(namef))
              difz = str('dz')
              dxdzf = str('{}/{}'.format(difnamef,difz))
              eqnf = str('{}{}{}{}{}{}'.format('{',difnamef,'}','{','dz','}'))
              d4 = Latex(r'$\frac{}$'.format(eqnf))
              d4b = Latex(r'${}$'.format(latexvalf))
              eqnvalf = str("{} = {}".format(sp.latex(sp.sympify(dxdzf)),latexvalf))
              latexeqnsf.append(eqnvalf)
       elif i2 == 24:
              indexkf = i2
              namef = namesf[i2]
              latexvalf = sp.latex(j2)
              difnamef = str('dT')
              difz = str('dz')
              dxdzf = str('{}/{}'.format(difnamef,difz))
              eqnf = str('{}{}{}{}{}{}'.format('{',difnamef,'}','{','dz','}'))
              d4 = Latex(r'$\frac{}$'.format(eqnf))
              d4b = Latex(r'${}$'.format(latexvalf))
              eqnvalf = str("{} = {}".format(sp.latex(sp.sympify(dxdzf)),latexvalf))
              latexeqnsf.append(eqnvalf)
       elif i2 == 25:
              indexkf = i2
              namef = namesf[i2]
              latexvalf = sp.latex(j2)
              difnamef = str('d^2 T')
              difz = str('d z^2')
              dxdzf = str('{}/{}'.format(difnamef,difz))
              dd2 = "{}{}{}{}{}{}{}".format("\frac",'{',difnamef,'}','{','dz','}')
              eqnf = str('{}{}{}{}{}{}'.format('{',difnamef,'}','{',difz,'}'))
              d4 = Latex(r'$\frac{}$'.format(eqnf))
              d4b = Latex(r'${}$'.format(latexvalf))
              eqnvalf = str(r'\frac{d^2 T}{dz^2} = '+ r'{}'.format(latexvalf))
              latexeqnsf.append(eqnvalf)
xvals = []
yvals = []
x,y = Jac.shape
for ix in range(x):
       for iy in range(y):
              val = Jac[ix,iy]
              val2 = sp.latex(val)
              yvals.append(val2)
       xvals.append(yvals[:])
       yvals.clear()
rcParams["text.usetex"] = False
fig = plt.figure(figsize=(x+1, y+1))
for krow, row in enumerate(xvals):
    for kcol, num in enumerate(row):
        plt.text(10*kcol + 15, 10*krow + 15, r'${}$'.format(num),   horizontalalignment='center',  verticalalignment='center',size = 5)
plt.axis([0, 10*(x + 1), 10*(y + 1), 0],'Image')
plt.xticks(np.linspace(0, 10*(x + 1), x + 2), [])
plt.yticks(np.linspace(0, 10*(y + 1), y + 2), [])
plt.grid(linestyle="solid")
#fig.tight_layout()
plt.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\Jacobian.svg', dpi=300)
plt.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\Jacobian.pdf', dpi=300)
d7 = sp.latex(Jac)
d4b = sp.Matrix([*ydot])
d5 = []
for i in ydot:
       val5 = sp.latex(i)
       d5.append(val5)
llen = len(d5)
fig = plt.figure(figsize=(llen,llen*2.0)) #figsize=(llen, llen)

intnum = 0
for i2 in d5:
#       ax = fig.add_axes([0,0,1,1])
       left, width = .25, .5
       bottom, height = .25, .5
       right = left + width
       top = bottom + height
#       ax.set_axis_off()

       listvl = llen - intnum
       indexh = d5.index(i2)
       nameval = str('{}'.format(nameslatex[indexh]))
       nameval2 = str('{}{}{}'.format('{',nameval,'}'))
       text = fig.text(0.0, 0.0075*listvl,r'$r_{}$ = '.format(nameval2) + r'${}$'.format(d5[indexh]),size = 10,wrap=False, multialignment='left')
       intnum += 1
fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\Ydot.svg', bbox_inches='tight')
fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\Ydot.pdf', bbox_inches='tight')
intnum2 = 0
for i2 in d5:
       fig = plt.figure()
       ax = fig.add_axes([0,0,1,1])
       left, width = .25, .5
       bottom, height = .25, .5
       right = left + width
       top = bottom + height
       ax.set_axis_off()
       listvl = llen - intnum2
       indexh = d5.index(i2)
       nameval = str('{}'.format(nameslatex[indexh]))
       nameval2 = str('{}{}{}'.format('{',nameval,'}'))
       text = ax.text(0.5*(left+right), 0.5*(bottom+top),r'$r_{}$ = '.format(nameval2) + r'${}$'.format(d5[indexh]), va= 'center',ha= 'center', bbox= dict(boxstyle="round", fc="white", alpha= 0.3, ec="black", pad=0.2))
       fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\r_{}.svg'.format(nameval), bbox_inches='tight')
       fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Reactions\r_{}.pdf'.format(nameval), bbox_inches='tight')
       intnum2 += 1

#print(Jac)
#ydot/ydot2 are the symbolic equations
#ydot2a/ydot2b are the symbolic equations with values for Ea and k_0 plugged in
#ys,ysf are symbolic representations of the chemical species 
#kd,kdf are symbolic representations of the K Values
#EA,EAf are symbolic representations of the Activation Energyies
#Jac/Jacf or the analytical jacobians, Jacf simply has more numbers plugged in
#ydot2,ydot2b, ysf, kdf, EAf, IMf, Jac2,Jac2f = makesymbsf(Initreactionsf,Eqlistf,namesf,tempsf,consf)  
#latexeqnsf = []
#for i in ydot2:
#       indexk = ydot2.index(i)
#       namef = namesf[indexk]
#       latexval = sp.latex(i)
#       difname = str('dC_{}'.format(namef))
#       difz = str('dz')
#       dxdz = str('{}/{}'.format(difname,difz))
#       eqn = str('{}{}{}{}{}{}'.format('{',difname,'}','{','dz','}'))
#       d4 = Latex(r'$\frac{}$'.format(eqn))
#       d4b = Latex(r'${}$'.format(latexval))
##       d5 = Latex(r'$\frac{}{}{}{}{}{}$ = '.format('{',difname,'}','{','dz','}') + r'${}$'.format(i))
#       eqnval = str("{} = {}".format(sp.latex(sp.sympify(dxdz)),latexval))
#       latexeqns.append(eqnval)
##       print(eqn)
##       print(r'$\frac{}{}{}{}{}{}$ = '.format('{',difname,'}','{','dz','}') + r'${}$'.format(i))
#print(latexeqnsf[0])
xvalsf = []
yvalsf = []
xf,yf = Jac2.shape
for ixf in range(xf):
       for iyf in range(yf):
              val = Jac2[ixf,iyf]
              val2 = sp.latex(val)
              yvalsf.append(val2)
       xvalsf.append(yvalsf[:])
       yvalsf.clear()
rcParams["text.usetex"] = False
fig = plt.figure(figsize=(xf+1, yf+1))
for krow, row in enumerate(xvalsf):
    for kcol, num in enumerate(row):
        plt.text(10*kcol + 15, 10*krow + 15, r'${}$'.format(num),   horizontalalignment='center',  verticalalignment='center',size = 5)
plt.axis([0, 10*(xf + 1), 10*(yf + 1), 0],'Image')
plt.xticks(np.linspace(0, 10*(xf + 1), xf + 2), [])
plt.yticks(np.linspace(0, 10*(yf + 1), yf + 2), [])
plt.grid(linestyle="solid")
#fig.tight_layout()
plt.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\FullJacobian.svg', dpi=300)
plt.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\FullJacobian.pdf', dpi=300)
#d7 = sp.latex(Jac2)
#d4bf = sp.Matrix([*ydot2])
d5f = []
for i in ydot2:
       val5 = sp.latex(i)
       d5f.append(val5)
llenf = len(d5f)
fig = plt.figure(figsize=(llenf,llenf*2.0)) #figsize=(llen, llen)

intnumf = 0
for i2 in d5f:
#       ax = fig.add_axes([0,0,1,1])
       left, width = .25, .5
       bottom, height = .25, .5
       right = left + width
       top = bottom + height
#       ax.set_axis_off()

       listvl = llenf - intnumf
       indexh = d5f.index(i2)
       nameval = str('{}'.format(nameslatexf[indexh]))
       nameval2 = str('{}{}{}'.format('{',nameval,'}'))
       text = fig.text(0.0, 0.0075*listvl,r'$r_{}$ = '.format(nameval2) + r'${}$'.format(d5f[indexh]),size = 10,wrap=False, multialignment='left')
       intnumf += 1
fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\Ydotfull.svg', bbox_inches='tight')
fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\Ydotfull.pdf', bbox_inches='tight')
intnum2f = 0
for i2 in d5f:
       fig = plt.figure()
       ax = fig.add_axes([0,0,1,1])
       left, width = .25, .5
       bottom, height = .25, .5
       right = left + width
       top = bottom + height
       ax.set_axis_off()
       listvl = llenf - intnum2f
       indexh = d5f.index(i2)
       nameval = str('{}'.format(nameslatexf[indexh]))
       nameval2 = str('{}{}{}'.format('{',nameval,'}'))
       text = ax.text(0.5*(left+right), 0.5*(bottom+top),r'$r_{}$ = '.format(nameval2) + r'${}$'.format(d5f[indexh]), va= 'center',ha= 'center', bbox= dict(boxstyle="round", fc="white", alpha= 0.3, ec="black", pad=0.2))
       fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\r_{}full.svg'.format(nameval), bbox_inches='tight')
       fig.savefig(r'C:\Users\tjcze\Desktop\Thesis\Python\Full\Reactions\r_{}full.pdf'.format(nameval), bbox_inches='tight')
       intnum2f += 1

def passarg(Jac,kd,EA,eav,kdv,cons,cons1,Twallval):
    k_0v = [*kdv]
    EAv = [*eav]
    allargs = k_0v
    for i in eav:
        allargs.append(i)
    Eas = []
    kv = []
    for i,ii in enumerate(EA):
        Eas.append(sp.symbols('{}'.format(EA[i])))
    for j,jj in enumerate(kd):
        kv.append(sp.symbols('{}'.format(kd[j])))
    kdict = {kd[k] : k_0v[k] for k in range(len(kd))}
    eadict = {Eas[l] : De(EAv[l]*1000) for l in range(len(EA))}
    rvald = {'Rval' : 8.314, 'Twall': Twallval, 'con1' : cons1[0], 'con2' : cons1[1], 'con3' : cons1[2]}
    dall = {}
    for i in (kdict,eadict,rvald):
        dall.update(kdict)
        dall.update(eadict)
        dall.update(rvald)
    Jaca = Jac.subs(dall)
    Jacd = sp.matrix2numpy(Jaca)
    Jace = sp.Matrix(Jacd)
    Jacf2 = Jace.xreplace({Rval: De(Rvalc)}) #
    Jf3 = sp.matrix2numpy(Jacf2)
    J4 = sp.Matrix(Jf3)
    return J4

def passargmp(Jac,kd,EA,eav,kdv,cons):
    Rval = sp.symbols('Rval')
    Twall = sp.symbols('Twall')
    k_0v = [*kdv]
    EAv = [*eav]
    allargs = k_0v
    con1 = consf[0]
    con2 = consf[1]
    con3 = consf[2]
    for i in eav:
        allargs.append(i)
    Eas = []
    kv = []
    for i,ii in enumerate(EA):
        Eas.append(sp.symbols('{}'.format(EA[i])))
    for j,jj in enumerate(kd):
        kv.append(sp.symbols('{}'.format(kd[j])))
    kdict = {kd[k] : k_0v[k] for k in range(len(kd))}
    eadict = {Eas[l] : mpf(EAv[l]*1000) for l in range(len(EA))}
    rvald = {Rval: mpf(8.314), Twall: mpf(Twallv), con1 : cons[0], con2 : cons[1], 'con3' : cons[2]}
    dall = {}
    for i in (kdict,eadict,rvald):
        dall.update(kdict)
        dall.update(eadict)
        dall.update(rvald)
    Jaca = Jac.subs(dall)
    Jacd = sp.matrix2numpy(Jaca)
    Jace = sp.Matrix(Jacd)
    Jacf2 = Jace.xreplace({Rval: mpf(Rvalc)}) #
    Jf3 = sp.matrix2numpy(Jacf2)
    J4 = sp.Matrix(Jf3)
    return J4

def jacobiana(symbolsjac,Jac):
    J_func = sp.lambdify(symbolsjac,Jac,modules='numpy')
    return J_func
    
def jacobianmp(symbolsjac,Jac):
    J_funcd = sp.lambdify(symbolsjac,Jac,modules='mpmath')
    return np.array([J_funcd],dtype=np.dtype(mpf))

args2 = [5.9E15, 1.3E13, 1.0E12, 5.0E11, 1.2E13, 2.0E11, 1.0E13, 1.0E13, 1.7E13, 1.2E13, 1.7E13, 9.1E10, 1.2E14, 5.0E11, 2.0E10, 3.0E11, 2.1E14, 5.0E14, 2.0E13, 1.0E14, 1.6E14] #If first order: [1/s]. If second order: [cm^3/mol*s]
args = [342.0,7.0,42.0,45.0,34.0,48.0,13.0,12.0,4.0,6.0,15.0,0.0,56.0,31.0,30.0,61.0,84.0,90.0,70.0,20.0,70.0]
argsf = [5.90E15,	2.20E12,	1.30E13,	1.20E13,	1.00E12,	5.00E11,	2.00E11,	1.00E11,	1.00E12,	1.00E13,	1.00E13,	1.70E13,	1.20E13,	1.70E13,	1.70E13,	1.60E13,	9.10E10,	1.20E14,	3.00E11,	2.00E10,	5.00E11,	2.10E14,	5.00E14,	2.00E13,	2.50E13,	1.00E12, 5.00E11,	5.00E11,	1.00E13,	1.00E14,	1.60E14]
args2f = [342,230,7,34,42,45,48,56,63,13,12,4,6,15,17,14,0,56,61,30,31,84,90,70,70,33,33,33,13,20,70]

#print(Jac)
#JacU = passarg(Jac,kd,EA,args,args2,consf,consv,Twallv)
#pp.pprint(JacU)

#JacUf = passarg(Jac2,kdf,EAf,argsf,args2f,cons)
#Jac_f = jacobiana(ys,JacU)
#JacUfmp = passargmp(Jac2,kdf,EAf,argsf,args2f,cons)
#Jac_fmp = jacobiana(ysf,JacUfmp)
#print(Jac_f)
#preview(JacUfmp,output='png',filename='test.png')
#Jf = np.asarray(JacUfmp)
#print(Jf)
#xx = len(ydot)
#new_img = Image.new("L", (xx, 1), "white")
#new_img.putdata(ydot)
#new_img.save(r'C:\Users\tjcze\Desktop\Thesis\Python\Simple\Equations2.png')
#img = Image.fromarray(ydot, 'RGB')
#img.save('Equations.png')
#new_img.show()


