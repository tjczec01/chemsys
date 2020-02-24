# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 09:11:20 2020

@author: tjcze
"""

import sympy as sp
import matplotlib.pyplot as plt 
from IPython.display import Latex

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
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" : 1, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C2H2',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :-1, "Reaction 19" :0, "Reaction 20" : -1, "Reaction 21" :-1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :-1, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C11',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : -1, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'C112',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 0, "Reaction 11" :-1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R1',
 "Reactions" : {"Reaction 1" : 1 ,"Reaction 2" : -1,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : -1, "Reaction 8" : -1, "Reaction 9" : -1, "Reaction 10" : -1, "Reaction 11" : -1, "Reaction 12" : -1, "Reaction 13" : -1, "Reaction 14" : 0, "Reaction 15" : 1, "Reaction 16" : 1, "Reaction 17" :1, "Reaction 18" :1, "Reaction 19" : 1, "Reaction 20" : 1, "Reaction 21" : -1}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :-1, "Reaction 18" :-1, "Reaction 19" :-1, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R2',
 "Reactions" : {"Reaction 1" : 1 ,"Reaction 2" : 0,"Reaction 3" : -1,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" : -1, "Reaction 8" :0, "Reaction 9" : 1, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : -1, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R3',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 1,"Reaction 3" : 1,"Reaction 4" : 1,"Reaction 5" : 1,"Reaction 6" : 1, "Reaction 7" :0, "Reaction 8" : -1, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" : -1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R4',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : -1,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" : 1, "Reaction 11" :0, "Reaction 12" : 1, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" : -1, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R5',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" : -1, "Reaction 14" : -1, "Reaction 15" :0, "Reaction 16" :1, "Reaction 17" :0, "Reaction 18" :-1, "Reaction 19" :0, "Reaction 20" :1, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :1, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'R6',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : -1, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" : 1, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" :0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :-1, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :0, "Reaction 18" :0, "Reaction 19" :1, "Reaction 20" :0, "Reaction 21" :0}},
{"Name" : 'VCM',
 "Reactions" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 1,"Reaction 6" : 0, "Reaction 7" :1, "Reaction 8" :0, "Reaction 9" :0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" : -1, "Reaction 13" : -1, "Reaction 14" : -1, "Reaction 15" : -1, "Reaction 16" : -1, "Reaction 17" : 1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}, 
 "Reverse" : {"Reaction 1" : 0 ,"Reaction 2" : 0,"Reaction 3" : 0,"Reaction 4" : 0,"Reaction 5" : 0,"Reaction 6" : 0, "Reaction 7" :0, "Reaction 8" :0, "Reaction 9" : 0, "Reaction 10" :0, "Reaction 11" :0, "Reaction 12" :0, "Reaction 13" :0, "Reaction 14" : 0, "Reaction 15" :0, "Reaction 16" :0, "Reaction 17" :-1, "Reaction 18" :0, "Reaction 19" :0, "Reaction 20" :0, "Reaction 21" :0}}, 
]


def symfunc(names, rxnum):
       Csyms = [sp.symbols(r'C_{}'.format('{}'.format(i))) for i in names]
       Ksyms = [sp.symbols(r'K_{}'.format(j)) for j in range(rxnum)]
       EAsyms = [sp.symbols(r'Ea_{}'.format(k)) for k in range(rxnum)]
       Tsyms = [sp.symbols('T')]
       
       return Csyms, Ksyms, EAsyms, Tsyms

def rterm(Ci, α):
       termi = sp.Mul(α,sp.Pow(Ci, abs(int(α))))
       return termi

def rprod(Ci, α, Cj, β):
       term1 = rterm(Ci, α)
       term2 = rterm(Cj, β)
       term3 = sp.Mul(term1,term2)
       return term3

def init(initlist, C):
       reactants = []
       products = []
       for i, j in enumerate(initlist):
              Reactants = initlist[i]['Reactants']
              Products = initlist[i]['Products']
              Rvals = list(Reactants.values())
              Pvals = list(Products.values())
              Ks = sp.symbols('K_{}'.format(i+1))
              Eas = sp.symbols('Ea_{}'.format(i+1))
              ee = sp.exp(sp.Mul(Eas,sp.Pow(sp.Mul(sp.Symbol('R'),sp.Symbol('T')), sp.Integer(-1))))
              rterms = []
              pterms = []
              rtotal = sp.Integer(1)
              ptotal = sp.Integer(1)
              for k, l in zip(C,Rvals):
                     if l != 0:
                            term = rterm(k, l)
                            rterms.append(term)
              for t in rterms:
                     rtotal = sp.Mul(rtotal,t)
              for m, n in zip(C,Pvals):
                     if n != 0:
                            pterm = rterm(m, n)
                            pterms.append(pterm)
              for tt in pterms:
                     ptotal = sp.Mul(ptotal,tt)
              reactants.append(sp.Mul(Ks,sp.Mul(rtotal,ee)))
              products.append(sp.Mul(Ks,sp.Mul(ptotal,ee)))
       return [reactants, products]

def eqlist(eqlistl, R, P):
       reactants = R
       products = P
       EQS = []
       leqns = []
       for i ,j in enumerate(eqlistl):
              Reactions = eqlistl[i]['Reactions']
              Reverse = eqlistl[i]['Reverse']
              Rxn = list(Reactions.values())
              RxnR = list(Reverse.values())
              eqn = []
              Reacts = [i*j for i, j in zip(Rxn, reactants) if i != 0]
              Prods = [i*j for i, j in zip(RxnR, products) if i != 0]
              if not Prods:
                     eee= sum(Reacts)
                     rlatex = sp.latex(eee)
                     leqns.append(rlatex)
                     EQS.append(eee)
              else:
                     eqn = sum(Reacts)
                     peqn = sum(Prods)
                     eeqn = sp.Add(eqn,peqn)
                     rlatex = sp.latex(eeqn)
                     leqns.append(rlatex)
                     EQS.append(eeqn)
       return [EQS, leqns]

def dislat(lnames, latexs, indvar):
       Latexs = []
       Displays = []
       for i in range(len(latexs)):
              dd = '{d'+ '{}'.format(sp.symbols(lnames[i]))+ '}'
              dt = '{d'+ '{}'.format(sp.symbols(indvar))+ '}'
              dde = r'$\frac{}{}'.format(dd,dt) + ' = ' + '{}$'.format(latexs[i])
              ddg = Latex(dde)
              Latexs.append(dde)
              Displays.append(ddg)
       return Displays, Latexs

def rhseqs(equations):
       EQLIST = []
       for ind, e in enumerate(equations):
              eqn = [r'{}'.format(e).replace('*exp',']*sp.exp').replace('{','').replace('}','').replace('K_','K[').replace('Ea_','Ea[').replace('/',']/')]
              for ind2 in range(len(equations)):
                     k1 = 'K[{}]'.format(ind2+1)
                     k2 = 'K[{}]'.format(ind2)
                     K1 = str(k1)
                     K2 = str(k2)
                     e1 = 'Ea[{}]'.format(ind2+1)
                     e2 = 'Ea[{}]'.format(ind2)
                     E1 = str(e1)
                     E2 = str(e2)
                     eqn2 =  str(eqn[0]).replace(E1,E2).replace(K1,K2)
                     eqn[0] = eqn2
              EQLIST.append(eqn[0])
       return EQLIST
          
indvar = str(input('Input independent variable --> '))    
number_of_reactions = len(Initreactions)
nameslist = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM']
Csymbols, K, EA, T = symfunc(nameslist, number_of_reactions)
reactants, products = init(Initreactions, Csymbols)
equations, latexs = eqlist(Eqlist, reactants, products)
slatex, dlatex = dislat(nameslist, latexs, indvar)
reqs = rhseqs(equations)
try:
       pvar = str(input('Print y/n --> '))
       if pvar == str('y'):
              printvar = True
       elif pvar == str('n'):
              printvar = False
       if printvar is True:
              filename = input('Enter save path --> ')
       else: 
              pass
except:
       print("Value Error: Input must be either 'y' or 'n'")
       

if printvar == True:
       for s,k in enumerate(dlatex):
              fig = plt.figure()
              ax = fig.add_axes([0,0,1,1])
              left, width = .25, .5
              bottom, height = .25, .5
              right = left + width
              top = bottom + height
              ax.set_axis_off()
              text = ax.text(0.5*(left+right), 0.5*(bottom+top),k , va= 'center',ha= 'center', bbox= dict(boxstyle="round", fc="white", alpha= 0.3, ec="black", pad=0.2))
              fig.savefig(r'{}\{}.svg'.format(filename,nameslist[s]), bbox_inches='tight')
              fig.savefig(r'{}\{}.pdf'.format(filename,nameslist[s]), bbox_inches='tight')
