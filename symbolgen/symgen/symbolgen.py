# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 01:14:06 2020

@author: tjcze
"""

import scipy as sc
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
# import symbolgen as sg
from tkinter import Tk, IntVar, StringVar, ttk, N, W, E, S, Checkbutton
# from tkinter import *

chemnumsl = []
rxnsvl = []
chemnamesl = []
reactants_num = []
products_num = []
reverse = []
coeffsr = []
coeffsp = []
Initreactions = []
Eqlist = []

def close_window():
     global entry
     entry = int(chems.get())
     chemnumsl.append(entry)
     entry2 = int(rxns.get())
     rxnsvl.append(entry2)
     root.destroy()

def close_window3():
     global entry
     for i in range(0,chemnumsl[0],1):
            entry3 = str(entries[i].get())
            chemnamesl.append(entry3)
     root3.destroy()

def close_window4():
     global entry
     for i in range(0,rxnsvl[0],1):
            entry4a = int(entriesr[i].get())
            reactants_num.append(entry4a)
            entry4b = int(entriesp[i].get())
            products_num.append(entry4b)
            entry4c = int(intvars[i].get())
            reverse.append(entry4c)
     root4.destroy()

def close_window5():
     global entry
     num_chems = int(len(chemnamesl))
     for i in range(0,rxnsvl[0],1):
            cfsr = [0*ij for ij in range(0,num_chems,1)]
            cfsp = [0*ik for ik in range(0,num_chems,1)]
            for j in range(0, reactants_num[i],1):
                   entry5r = entriesrc[i][j].get()
                   indexr = chemnamesl.index(entry5r)
                   cfsr[indexr] = int(entriesr5[i][j].get())
            coeffsr.append(cfsr[:]) 
            cfsr.clear()
            for k in range(0, products_num[i],1):
                   entry5p = entriespc[i][k].get()
                   indexp = chemnamesl.index(entry5p)
                   cfsp[indexp] = int(entriesp5[i][k].get())
            coeffsp.append(cfsp[:]) 
            cfsp.clear()
     root5.destroy()

root = Tk()
root.title("Number of chemical species")
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)	
chemnums = StringVar()
chems = ttk.Entry(mainframe, width=7, textvariable=chemnums)
chems.grid(column=2, row=1, sticky=(W, E))
ttk.Label(mainframe, text="Enter total number of chemical species ").grid(column=1, row=1, sticky=W)
rnums = StringVar()
rxns = ttk.Entry(mainframe, width=7, textvariable=rnums)
rxns.grid(column=2, row=2, sticky=(W, E))
ttk.Label(mainframe, text="Enter total number of chemical reactions ").grid(column=1, row=2, sticky=W)
B = ttk.Button(root, text = "OK", command = close_window).grid(column=3, row=1)
root.mainloop()

chems_value = chemnumsl[0]

rxnnum = rxnsvl[0]

root3 = Tk()
root3.title("Name of chemical species")
mainframe3 = ttk.Frame(root3, padding="3 3 12 12")
mainframe3.grid(column=0, row=0, sticky=(N, W, E, S))
root3.columnconfigure(0, weight=1)
root3.rowconfigure(0, weight=1)	
stringvars = []
entries = []
for i in range(0,chems_value,1):
       stringvars.append(StringVar())
for i in range(0,chems_value,1):
       entries.append(ttk.Entry(mainframe3, width=20, textvariable=stringvars[i]))
       entries[i].grid(column=2, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe3, text="Enter name of chemical species {} ".format(i+1)).grid(column=1, row=int(i+1), sticky=W)
B3 = ttk.Button(root3, text = "OK", command = close_window3).grid(column=3, row=1)
root3.mainloop()

root4 = Tk()
root4.title("Number of reactants & products")
mainframe4 = ttk.Frame(root4, padding="3 3 12 12")
mainframe4.grid(column=0, row=0, sticky=(N, W, E, S))
root4.columnconfigure(0, weight=1)
root4.rowconfigure(0, weight=1)	
stringvarsr = []
stringvarsp = []
intvars = []
entriesr = []
entriesp = []
entriesc = []
for i in range(0,rxnnum,1):
       stringvarsr.append(StringVar())
       stringvarsp.append(StringVar())
       intvars.append(IntVar())
for i in range(0,rxnnum,1):
       entriesr.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsr[i]))
       entriesr[i].grid(column=2, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text=" Enter number of reactants in reaction {} ".format(i+1)).grid(column=1, row=int(i+1), sticky=W)
       entriesp.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsp[i]))
       entriesp[i].grid(column=4, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text=" Enter number of products in reaction {} ".format(i+1)).grid(column=3, row=int(i+1), sticky=W)
       entriesc.append(Checkbutton(root4, text=" Reversable", variable=intvars[i]))
       entriesc[i].grid(column=5, row=int(i+1), sticky=W)
B4 = ttk.Button(root4, text = "OK", command = close_window4).grid(column=2, row=int(rxnnum+2))
root4.mainloop()

root5 = Tk()
root5.title("Reactions")
mainframe5 = ttk.Frame(root5, padding="3 3 12 12")
mainframe5.grid(column=0, row=0, sticky=(N, W, E, S))
root5.columnconfigure(0, weight=1)
root5.rowconfigure(0, weight=1)	
stringvarsr5 = []
stringvarsp5 = []
stringvarsrc = []
stringvarspc = []
entriesr5a = []
entriesp5a = []
entriesr5 = []
entriesp5 = []
entriesrca = []
entriespca = []
entriesrc = []
entriespc = []
rstrings = []

num_chems = len(chemnamesl) 
for i in range(0,rxnnum,1):
       rval = reverse[i]
       if rval == 0:
              rstrings.append(str(' --> '))
       elif rval == 1:
              rstrings.append(str(' <==> '))
for i in range(0,rxnnum,1):
       stringvarsr5.append([StringVar() for i in range(0,reactants_num[i],1)])
       stringvarsp5.append([StringVar() for i in range(0,products_num[i],1)])
       stringvarsrc.append([StringVar() for i in range(0,reactants_num[i],1)])
       stringvarspc.append([StringVar() for i in range(0,products_num[i],1)])

for i in range(0,rxnnum,1):
       rangeval = int(reactants_num[i])
       rangep = int(products_num[i])
       rval2 = int(reactants_num[i])*2 + 1
       int1 = 1
       int2 = 2
       jval = 1
       for j in range(0,reactants_num[i],1):
              entriesr5a.append(ttk.Entry(mainframe5, width=7, textvariable=stringvarsr5[i][j]))
              entriesr5a[-1].grid(column=jval, row=int(i+1))
              jval += 1
              combbo = ttk.Combobox(mainframe5, values=chemnamesl)
              combbo.grid(column=jval, row=int(i+1))
              jval += 1
              entriesrca.append(combbo) 
              if j < reactants_num[i]-1:
                     ttk.Label(mainframe5, text=" + ").grid(column=jval, row=int(i+1))
                     jval += 1
              elif j == reactants_num[i]-1: 
                     ttk.Label(mainframe5, text=" {} ".format(rstrings[i])).grid(column=jval, row=int(i+1))
                     jval += 1
              int1 += 1
              int2 += 1
       endp = rval2 + products_num[i]+1
       entriesr5.append(entriesr5a[:])
       entriesr5a.clear()
       entriesrc.append(entriesrca[:])
       entriesrca.clear()
       for k in range(0, products_num[i],1):
              int1b = 1
              int2b = 2
              entriesp5a.append(ttk.Entry(mainframe5, width=7, textvariable=stringvarsp5[i][k]))
              entriesp5a[-1].grid(column=jval, row=int(i+1))
              jval += 1
              combbb = ttk.Combobox(mainframe5, values=chemnamesl)
              combbb.grid(column=jval, row=int(i+1))
              jval += 1
              entriespca.append(combbb)
              if k < products_num[i]-1:
                     ttk.Label(mainframe5, text=" + ").grid(column=jval, row=int(i+1))
                     jval += 1
              else: 
                     passB5 = ttk.Button(root5, text = "OK", command = close_window5).grid(column=2, row=int(rxnnum+2))
       entriesp5.append(entriesp5a[:])
       entriesp5a.clear()
       entriespc.append(entriespca[:])
       entriespca.clear()
root5.mainloop()

rxns_strs = ["Reaction {}".format(int(i+1)) for i in range(0,rxnnum,1)]

for i in range(0,rxnnum,1):
       indexnum = int(i+1)
       keys = chemnamesl
       valuesr = coeffsr[i][:]
       valuesp = coeffsp[i][:]
       dictionary = {"Ea" : indexnum, "K_Value": indexnum, "Reactants" : dict(zip(keys , valuesr)), "Products" : dict(zip(keys , valuesp))}
       Initreactions.append(dictionary)
       
for i in range(0,len(chemnamesl),1):
       indexnum = int(i+1)
       namev = chemnamesl[i]
       name_index = chemnamesl.index(namev)
       keys = rxns_strs
       valuesfor = [0*ij for ij in range(0,rxnnum,1)]
       valuesrev = [0*ik for ik in range(0,rxnnum,1)]
       for j in range(0,rxnnum,1):
              valuef = coeffsr[j][name_index]
              if valuef != 0:
                     valuesfor[j] = int(-1)
              elif coeffsp[j][name_index] != 0:
                     valuesfor[j] = int(1)
              valuesrev[j] = int(1*reverse[j]) + int(-1*reverse[j]) #coeffsp[j][name_index]*
       dictionary2 = {"Name" : "{}".format(str(namev)), "Reactions" : dict(zip(keys , valuesfor)), "Reverse" : dict(zip(keys , valuesrev))}
       Eqlist.append(dictionary2)
# ]

class symbolgen:
       
       indvar = str(input('Input independent variable --> '))
       
       def __init__(self, nlist,  Initlist, EQlist):
              self.nameslist = nlist
              self.rxnnum = len(self.nameslist)
              self.initlist = Initlist
              self.Eqlist = EQlist
       
       def initl(self):
               return self.initlist
        
       def latexin(self):
              latexs = self.eqlist(Eqlist, self.reactants, self.products)[1]
              return latexs
       
       def symsinit(self):
               return self.symfunc(self.nameslist,self.rxnnum)[0]
       
       def rinit(self):
               return self.initfunc(self.initreactions, self.C)[0]
       
       def pinit(self):
               return self.initfunc(self.initreactions, self.C)[1]
       
       initreactions = property(initl)
       C = property(symsinit, symsinit)
       reactants = property(rinit)
       products = property(pinit)
       latexs = property(latexin)
        
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
              term1 = symbolgen.rterm(Ci, α)
              term2 = symbolgen.rterm(Cj, β)
              term3 = sp.Mul(term1,term2)
              return term3
       
       def initfunc(initlist, C):
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
                                   term = symbolgen.rterm(k, l)
                                   rterms.append(term)
                     for t in rterms:
                            rtotal = sp.Mul(rtotal,t)
                     for m, n in zip(C,Pvals):
                            if n != 0:
                                   pterm = symbolgen.rterm(m, n)
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
       
       def sysgen(self):
              equations, latexs = self.eqlist(self.Eqlist, self.reactants, self.products)
              return equations
       
       def sysdis(self):
              equations, latexs = self.eqlist(self.Eqlist, self.reactants, self.products)
              slatex, dlatex = self.dislat(self.nameslist, self.latexs, self.indvar)
              return dlatex
       
       def dis(self):
              slatex, dlatex = self.dislat(self.nameslist, self.latexs, self.indvar)
              for i in slatex:
                     display(i)
                     
       def gen(names, inits, eqs):
              rxn = len(names)
              Cs, Ks, EAs, Ts = symbolgen.symfunc(names, rxn)
              reacts, prods = symbolgen.initfunc(inits, Cs)
              equats, latexss = symbolgen.eqlist(eqs, reacts, prods)
              slat, dlat = symbolgen.dislat(names, latexss, symbolgen.indvar)
              return Cs, equats, slat, dlat
       
       def psave(nameslist, LATEXD):
              filename = input('Enter save path --> ')
              for s,k in enumerate(LATEXD):
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

names = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM']

C, equations, slat, dlat = sg.gen(chemnamesl, Initreactions, Eqlist)
sg.psave(names, slat)

