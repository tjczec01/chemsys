import sympy as sp
import matplotlib.pyplot as plt
from IPython.display import display, Latex
from tkinter import Tk, IntVar, StringVar, ttk, N, W, E, S, Checkbutton
import pickle

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
indvdf = []
ffpath = []
kk = []
eaf = []
RR = []
def close_window():
     global entry
     entry = int(chems.get())
     chemnumsl.append(entry)
     entry2 = int(rxns.get())
     rxnsvl.append(entry2)
     entry3 = str(indvard.get())
     indvdf.append(r'{}'.format(entry3))
     entry4 = str(r'{}'.format(filev.get()))
     ffpath.append(entry4)
     rval = float(rg.get())
     RR.append(rval)
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
            entryk = float(entriesk[i].get())
            kk.append(entryk)
            entryea = float(entriesea[i].get())
            eaf.append(entryea)
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
indvard = StringVar()
inv = ttk.Entry(mainframe, width=7, textvariable=indvard)
inv.grid(column=2, row=3, sticky=(W, E))
ttk.Label(mainframe, text="Enter independent variable ").grid(column=1, row=3, sticky=W)
filep = StringVar()
filev = ttk.Entry(mainframe, width=50, textvariable=filep)
filev.grid(column=2, row=4, sticky=(W, E))
ttk.Label(mainframe, text="Enter file path ").grid(column=1, row=4, sticky=W)
rgas = StringVar()
rg = ttk.Entry(mainframe, width=7, textvariable=rgas)
rg.grid(column=2, row=5, sticky=(W, E))
ttk.Label(mainframe, text="Enter Gas Constant ").grid(column=1, row=5, sticky=W)
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
stringvarsk = []
stringvarsea = []
intvars = []
entriesr = []
entriesp = []
entriesc = []
entriesk = []
entriesea = []
for i in range(0,rxnnum,1):
       stringvarsr.append(StringVar())
       stringvarsp.append(StringVar())
       stringvarsk.append(StringVar())
       stringvarsea.append(StringVar())
       intvars.append(IntVar())
for i in range(0,rxnnum,1):
       entriesr.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsr[i]))
       entriesr[i].grid(column=3, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text=" Enter number of reactants in reaction {} ".format(i+1)).grid(column=2, row=int(i+1), sticky=W)
       entriesp.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsp[i]))
       entriesp[i].grid(column=5, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text=" Enter number of products in reaction {} ".format(i+1)).grid(column=4, row=int(i+1), sticky=W)
       entriesk.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsk[i]))
       entriesk[i].grid(column=7, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text="k_{0}  = ".format(i+1)).grid(column=6, row=int(i+1), sticky=W)
       entriesea.append(ttk.Entry(mainframe4, width=7, textvariable=stringvarsea[i]))
       entriesea[i].grid(column=9, row=int(i+1), sticky=(W, E))
       ttk.Label(mainframe4, text="Ea_{0} [kJ/mol] = ".format(i+1)).grid(column=8, row=int(i+1), sticky=W)
       entriesc.append(Checkbutton(root4, text="Reaction {} Reversable".format(i+1), variable=intvars[i]).grid(column=10, row=int(i), sticky=(N,E)))
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
       dictionary = {"Ea" : indexnum, "K_Value": indexnum, "Reverse" : reverse[i], "Reactants" : dict(zip(keys , valuesr)), "Products" : dict(zip(keys , valuesp))}
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



class symbolgen:
       
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
       
       def numfunc(Cs):
              cl = len(Cs)
              fcs = []
              for i in range(cl):
                     As = []
                     Ns = []
                     NNs = []
                     val3 = Cs[i]
                     se = list(val3)
                     count = 0
                     sb = list(val3)
                     SG = []
                     fnum = len(se) - 1
                     fend = len(se) - 1
                     for sv in range(len(sb)):
                            fend = len(se)
                            vv = sb[sv]
                            N = vv.isnumeric()
                            A = vv.isalpha()
                            ff = fend - sv
                            if A == True and count == 0:
                                   As.append(vv)
                                   SG.append(vv)
                                   count = 0
                                   fnum -= 1
                            elif A == True and count > 0:
                                   
                                   NNa = "".join(Ns)
                                   SG.append(NNa)
                                   SG.append(vv)
                                   Ns.clear()
                                   count = 0
                                   fnum -= 1
                            elif A == True and count >= 2:
                                  
                                   NNa = "".join(Ns)
                                   NNs.append(NNa)
                                   Ns.clear()
                                   SG.append(NNa)
                                   SG.append(vv)
                                   count = 0
                                   fnum -= 1
                                   
                            elif N == True and ff > 1:
                                   Ns.append(vv)
                                   count += 1
                            elif N == True and ff <= 1:
                                   Ns.append(vv)
                                   
                                   if len(Ns) >= 2:
                                          NNa = "".join(Ns)
                                          NNs.append(NNa)
                                          SG.append(NNa)
                                                 
                                   else:
                                          SG.append(vv)
                                   
                     count = 0
                     Ns.clear()
                     As.clear()
                     val2 = str(Cs[i])
                     s = list(val2)
                     for j in range(len(SG)):
                            charv = SG[j]
                            try:
                                    charvi = int(SG[j])
                                    SG[j] = charv.replace('{}'.format(charvi), ('_{' + '{}'.format(charvi) + '}')) 
                            except:
                                    pass                                                 
                     ss = "".join(SG)
                     s.clear()
                     fcs.append(ss)
              return fcs
                            
                     
       
       def rterm(Ci, a):
              termi = sp.Mul(a,sp.Pow(Ci, abs(int(a))))
              return termi
       
       def rprod(Ci, a, Cj, b):
              term1 = symbolgen.rterm(Ci, a)
              term2 = symbolgen.rterm(Cj, b)
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
              return reactants, products
       
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
              return EQS, leqns
       
       def dislat(lnames, latexs, indvar):
              Latexs = []
              Displays = []
              for i in range(len(latexs)):
                     dd = '{d'+ '{}'.format(sp.symbols(lnames[i]))+ '}'
                     dt = '{d'+ '{}'.format(sp.symbols(indvar))+ '}'
                     dde = r'$\dfrac{}{}'.format(dd,dt) + ' = ' + '{}$'.format(latexs[i])
                     ddg = Latex(dde)
                     Latexs.append(dde)
                     Displays.append(ddg)
              return Displays, Latexs
       
       def chemeq(Cs, rxn, inits):
              ceqs = []
              ceqsD = []
              ceqsw = []
              for i in range(rxn):
                     Reactants = inits[i]['Reactants']
                     Products = inits[i]['Products']
                     Reverse = inits[i]['Reverse']
                     Rvals = list(Reactants.values())
                     rvals = [Rvals[kk] for kk in range(len(Rvals)) if Rvals[kk] != 0]
                     Rname = symbolgen.numfunc(list(Reactants.keys()))
                     rname = [sp.symbols('{}'.format(Rname[h])) for h in range(len(Rname)) if Rvals[h] !=0]
                     Pvals = list(Products.values())
                     pvals = [Pvals[kk] for kk in range(len(Pvals)) if Pvals[kk] != 0]
                     Pname = symbolgen.numfunc(list(Products.keys()))
                     pname = [sp.symbols('{}'.format(Pname[h])) for h in range(len(Pname)) if Pvals[h] !=0]
                     CRvals = sum([sp.Mul(sp.Integer(ii),jj) for ii, jj in zip(rvals,rname) if ii != 0])
                     CPvals = sum([sp.Mul(sp.Integer(ii),jj) for ii, jj in zip(pvals,pname) if ii != 0])
                     if Reverse == 0:
                            cheme = r'${} \longrightarrow {}$'.format(CRvals, CPvals) 
                     if Reverse == 1:
                            cheme = r'${} \rightleftharpoons {}$'.format(CRvals, CPvals) 
                     ceqsD.append(Latex(cheme))
                     ceqs.append(cheme)
                     if Reverse == 0:
                            chemw = r'{} \\longrightarrow {}'.format(CRvals, CPvals) 
                     if Reverse == 1:
                            chemw = r'{} \\rightleftharpoons {}'.format(CRvals, CPvals)
                     ceqsw.append(chemw)
              return ceqs, ceqsD, ceqsw
       
       def rhseqs(equations, kk, ea, r):
              EQLIST = []
              EQLISTF = []
              for ind, e in enumerate(equations):
                     eqn = [r'{}'.format(e).replace('{','').replace('}','')] 
                     Ksyms = [sp.symbols('K_{}'.format(i+1)) for i in range(len(kk))]
                     EAsyms = [sp.symbols('Ea_{}'.format(i+1)) for i in range(len(ea))]
                     kdictionary = dict(zip(Ksyms, kk))
                     eadictionary = dict(zip(EAsyms, ea))
                     eqn3 = e.subs(kdictionary)
                     eqn4 = eqn3.subs(eadictionary)
                     eqn5 = eqn4.subs({'R' : 8.31446261815324})
                     eqn6 = eqn5.subs({'*exp' : '*sp.exp'})
                     EQLISTF.append(eqn6)
                     EQLIST.append(eqn[0])
                     
              return EQLIST, EQLISTF
       
       def jacobian(rhs, y):
              funcs = [i for i in rhs]
              eqnl = len(funcs)
              cl = len(y)
              J  = [[i for i in range(cl)] for j in range(eqnl)]
              mfunc = lambda i, j: sp.diff(funcs[i], y[j])
              Jjj = sp.Matrix(eqnl, cl, lambda i, j: mfunc(i, j))
              Jj = sp.matrix2numpy(Jjj)
              for i in range(eqnl):
                      for j in range(cl):
                             J[i][j] = str('{}'.format('{}'.format(mfunc(i, j)).replace('*exp',  '*sp.exp')))
              MatrixJ = sp.Matrix(Jj)
              LatexMatrix = sp.latex(sp.matrix2numpy(sp.Matrix(Jj)))
              lm = sp.latex(MatrixJ, mat_delim="(" )
              return J, MatrixJ, lm, LatexMatrix
       
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
                     
       def gen(names, rxn, inits, eqs, intz):
              Cs, Ks, EAs, Ts = symbolgen.symfunc(names, rxn)
              reacts, prods = symbolgen.initfunc(inits, Cs)
              equats, latexss = symbolgen.eqlist(eqs, reacts, prods)
              slat, dlat = symbolgen.dislat(names, latexss, intz)
              Chem, ChemD, ChemW = symbolgen.chemeq(Cs, rxn, inits)
              return Cs, reacts, prods, equats, slat, dlat, Chem, ChemD, ChemW
       
       def fullgen(names, rxn, inits, eqs, intz, filepathf, kk, ea, r):
              Cs, Ks, EAs, Ts = symbolgen.symfunc(names, rxn)
              reacts, prods = symbolgen.initfunc(inits, Cs)
              equats, latexss = symbolgen.eqlist(eqs, reacts, prods)
              slat, dlat = symbolgen.dislat(names, latexss, intz)
              Chem, ChemD, ChemW = symbolgen.chemeq(Cs, rxn, inits)
              RHS, RHSf = symbolgen.rhseqs(equats, kk, ea, r)
              Jac, Jacm, lm, latexmatrix = symbolgen.jacobian(RHSf, Cs)
              fp = symbolgen.csave(Chem, filepathf)
              fc = symbolgen.psave(names, dlat, filepathf)
              ff = symbolgen.fsave(filepathf, equats, dlat, Chem, ChemW, RHS, RHSf, Jac, Jacm, lm, latexmatrix)
              return Cs, reacts, prods, equats, slat, dlat, Chem, ChemD, ChemW, RHS, RHSf, Jac, Jacm, lm, latexmatrix
       
       
       def psave(nameslist, LATEXD, fpath):
              filename = fpath
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
                     
       def csave(LATEXC, fpath):
              filename = fpath
              for s,k in enumerate(LATEXC):
                     fig = plt.figure()
                     ax = fig.add_axes([0,0,1,1])
                     left, width = .25, .5
                     bottom, height = .25, .5
                     right = left + width
                     top = bottom + height
                     ax.set_axis_off()
                     text = ax.text(0.5*(left+right), 0.5*(bottom+top),k , va= 'center',ha= 'center', bbox= dict(boxstyle="round", fc="white", alpha= 0.3, ec="black", pad=0.2))
                     fig.savefig(r'{}\Reaction {}.svg'.format(filename, s), bbox_inches='tight')
                     fig.savefig(r'{}\Reaction {}.pdf'.format(filename, s), bbox_inches='tight')
                     
       def fsave(ffpath, eqns, eqnslat, crxns, crxnsw, rhseq, rhseqf, Jac, JacM, lm, latexmatrix):
              with open(r"{}\Equations.txt".format(ffpath), "w") as output:
                     output.write("[")
                     for eqn in eqns:
                            output.write('{}'.format(str(eqn)))
                            output.write(",\n")
                     output.seek(0, 2)
                     output.seek(output.tell() - 2, 0)
                     output.write("]")
                  
              with open(r"{}\EquationsLatex.txt".format(ffpath), "w") as output:
                     for eqnlat in eqnslat:
                            output.write('{}\n'.format(str(eqnlat)))
                  
              with open(r"{}\ReactionsLatex.txt".format(ffpath), "w") as output:
                     for crxn in crxns:
                            output.write('{}\n'.format(str(crxn)))
                            
              with open(r"{}\ReactionsLatexWord.txt".format(ffpath), "w") as output:
                     for crxnw in crxnsw:
                            output.write('{}\n'.format(str(crxnw)))
                            
              with open(r"{}\RHSsymbols.txt".format(ffpath), "w") as output:
                     output.write("[")
                     for rhs in rhseq:
                            output.write('{}'.format(rhs))
                            output.write(",\n")
                     output.seek(0, 2)
                     output.seek(output.tell() - 2, 0)
                     output.write("]")
              
              with open(r"{}\RHS.txt".format(ffpath), "w") as output:
                     output.write("[")
                     for rhs in rhseqf:
                            output.write('{}'.format(rhs))
                            output.write(",\n")
                     output.seek(0, 2)
                     output.seek(output.tell() - 2, 0)
                     output.write("]")
                            
              with open(r"{}\Jacobun.txt".format(ffpath), "w") as output:
                     output.write("[")
                     for i in range(len(Jac)):
                             Jrow = Jac[i][:]
                             output.write('{}\n'.format(Jrow))
                     output.seek(0, 2)
                     output.seek(output.tell() - 2, 0)
                     output.write("]")
              
              with open("{}\JacobianMatrix.txt".format(ffpath),'w') as output:
                     output.write('{}'.format(JacM))
                            
              with open(r"{}\JacobianLatex.txt".format(ffpath), "w") as output:
                            output.write('{}'.format(lm))
                            
              with open("{}\Jacobun.txt".format(ffpath)) as filein, open("{}\Jacobian.txt".format(ffpath),'w') as fileout:
                     for line in filein:
                         line=line.replace("'","")
                         fileout.write(line)
                         
              with open(r"{}\RHS.txt".format(ffpath)) as filein, open(r"{}\RightHandSide.txt".format(ffpath),'w') as fileout:
                     for line in filein:
                         line=line.replace("'","")
                         line=line.replace("exp","sp.exp")
                         fileout.write(line)
                         
              with open(r"{}\RHSsymbols.txt".format(ffpath)) as filein, open(r"{}\RightHandSideSymbols.txt".format(ffpath),'w') as fileout:
                     for line in filein:
                         line=line.replace("'","")
                         line=line.replace("exp","sp.exp")
                         fileout.write(line)
                     
              pickle.dumps(JacM)      
              with open('{}\JacobianMatrixPickle.txt'.format(ffpath),'wb') as f:
                     pickle.dump(JacM, f)

names = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM']
ea = [i*1000.0 for i in eaf]
RRv = RR[0]
rxnnumf = rxnsvl[0]
C, reacts, prods, equations, slat, dlat, chem, chemD, chemw, rhs, rhsf, jac, jacM, lm, latexmatrix = symbolgen.fullgen(names, rxnnumf, Initreactions, Eqlist, indvdf[0], ffpath[0], kk, ea, RRv)
