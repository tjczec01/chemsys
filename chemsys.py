# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 01:14:06 2020

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

"""
import sympy as sp
from sympy import diff, Matrix, symbols, Add, Mul, Pow, Symbol, Integer, latex, exp, simplify
from sympy.matrices.dense import matrix2numpy
import matplotlib.pyplot as plt
from IPython.display import display, Latex
from tkinter import Tk, ttk, IntVar, StringVar, N, W, E, S, Checkbutton, Label, Entry, Button
from tkinter.ttk import Combobox
import pickle
import os
import subprocess
from shutil import which

def create_pdf(file_in, file_out):
    # print(which("pdflatex").replace("EXE", "exe") + ' -output-format=pdf ' + r"-output-directory={} ".format(file_out) + "-enable-pipes " + "-enable-mltex " + r"{}".format(file_in))
    cmds = str('"{}"'.format(which("pdflatex").replace("EXE", "exe") + ' -output-format=pdf ' + r"-output-directory={} ".format(file_out) + "-enable-pipes " + "-enable-mltex " + r"{}".format(file_in)))
    os.system(cmds)
    process = subprocess.Popen([which("pdflatex").replace("EXE", "exe"), '-output-format=pdf' , r"-output-directory={}".format(file_out) , "-enable-pipes" , "-enable-mltex" , r"{}".format(file_in)])
    process.wait()


clear = os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = r"{}\Jacobian".format(dir_path)

try:
    os.mkdir(path_fol)
except Exception:
    pass

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


def close_window2():
    global entry
    for i in range(0, chemnumsl[0], 1):
        entry2 = str(entries[i].get())
        chemnamesl.append(entry2)
    root2.destroy()


def close_window3():
    global entry
    for i in range(0, rxnsvl[0], 1):
        entry3a = int(entriesr[i].get())
        reactants_num.append(entry3a)
        entry3b = int(entriesp[i].get())
        products_num.append(entry3b)
        entry3c = int(intvars[i].get())
        entryk = float(entriesk[i].get())
        kk.append(entryk)
        entryea = float(entriesea[i].get())
        eaf.append(entryea)
        reverse.append(entry3c)
    root3.destroy()


def close_window4():
    global entry
    num_chems = int(len(chemnamesl))
    for i in range(0, rxnsvl[0], 1):
        cfsr = [0*ij for ij in range(0, num_chems, 1)]
        cfsp = [0*ik for ik in range(0, num_chems, 1)]
        for j in range(0, reactants_num[i], 1):
            entry4r = entriesrc[i][j].get()
            indexr = chemnamesl.index(entry4r)
            cfsr[indexr] = int(entriesr4[i][j].get())
        coeffsr.append(cfsr[:])
        cfsr.clear()
        for k in range(0, products_num[i], 1):
            entry4p = entriespc[i][k].get()
            indexp = chemnamesl.index(entry4p)
            cfsp[indexp] = int(entriesp4[i][k].get())
        coeffsp.append(cfsp[:])
        cfsp.clear()
    root4.destroy()


root = Tk()
root.title("Number of chemical species")
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
chemnums = StringVar()
chems = Entry(mainframe, width=7, textvariable=chemnums)
chems.grid(column=2, row=1, sticky=(W, N, E, S))
Label(mainframe, text="Enter total number of chemical species ").grid(column=1, row=1, sticky=(W, N, E, S))
rnums = StringVar()
rxns = Entry(mainframe, width=7, textvariable=rnums)
rxns.grid(column=2, row=2, sticky=(W, N, E, S))
Label(mainframe, text="Enter total number of chemical reactions ").grid(column=1, row=2, sticky=(W, N, E, S))
indvard = StringVar()
inv = Entry(mainframe, width=7, textvariable=indvard)
inv.grid(column=2, row=3, sticky=(W, N, E, S))
Label(mainframe, text="Enter independent variable ").grid(column=1, row=3, sticky=(W, N, E, S))
filep = StringVar(value=path_fol)
filev = Entry(mainframe, width=50, textvariable=filep)
filev.grid(column=2, row=4, sticky=(W, N, E, S))
Label(mainframe, text="Enter file path ").grid(column=1, row=4, sticky=(W, N, E, S))
rgas = StringVar(value="8.31446261815324")
rg = Entry(mainframe, width=7, textvariable=rgas)
rg.grid(column=2, row=5, sticky=(W, N, E, S))
Label(mainframe, text="Enter Gas Constant ").grid(column=1, row=5, sticky=(W, N, E, S))
B = Button(root, text="OK", command=close_window).grid(column=3, row=1)
root.mainloop()

chems_value = chemnumsl[0]

rxnnum = int(rxnsvl[0])
rxnnumr = [int(i + 1) for i in range(rxnnum)]

root2 = Tk()
root2.title("Name of chemical species")
mainframe2 = ttk.Frame(root2, padding="3 3 12 12")
mainframe2.grid(column=0, row=0, sticky=(N, W, E, S))
root2.columnconfigure(0, weight=1)
root2.rowconfigure(0, weight=1)
stringvars = []
entries = []
for i in range(0, chems_value, 1):
    stringvars.append(StringVar())
for i in range(0, chems_value, 1):
    entries.append(Entry(mainframe2, width=20, textvariable=stringvars[i]))
    entries[i].grid(column=2, row=int(i + 1), sticky=(W, N, E, S))
    Label(mainframe2, text="Enter name of chemical species {} ".format(i + 1)).grid(column=1, row=int(i + 1), sticky=(W, N, E, S))
B3 = Button(root2, text="OK", command=close_window2).grid(column=3, row=1)
root2.mainloop()

root3 = Tk()
root3.title("Reactants & Products")
mainframe4 = ttk.Frame(root3, padding="3 3 12 12")
mainframe4.grid(column=0, row=0, sticky=(N, W, E))
root3.columnconfigure(0, weight=1)
root3.rowconfigure(0, weight=1)
root3.rowconfigure(1, weight=1)
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
for i in rxnnumr:
    stringvarsr.append(StringVar())
    stringvarsp.append(StringVar())
    stringvarsk.append(StringVar())
    stringvarsea.append(StringVar())
    intvars.append(IntVar())
for i in rxnnumr:
    mainframe4.rowconfigure(i, weight=1)
    coli0 = 0
    coli1 = coli0 + 1
    coli2 = coli1 + 1
    coli3 = coli2 + 1
    coli4 = coli3 + 1
    coli5 = coli4 + 1
    coli6 = coli5 + 1
    coli7 = coli6 + 1
    coli8 = coli7 + 1
    coli9 = coli8 + 1
    coli10 = coli9 + 1
    clist = [i for i in range(coli10 + 1)]
    for ci in clist:
        mainframe4.columnconfigure(ci, weight=1)
    Pad_x = 5
    Pad_y = 2
    CE = 2
    group0 = Label(mainframe4, text="Reaction Number", padx=Pad_x, pady=Pad_y).grid(row=0, column=coli0, columnspan=1, sticky=(W, N, E, S))
    group1 = Label(mainframe4, text="Number of Reactants", padx=Pad_x, pady=Pad_y).grid(row=0, column=coli1, columnspan=2, sticky=(W, N, E, S))
    group2 = Label(mainframe4, text="Number of Products", padx=Pad_x, pady=Pad_y).grid(row=0, column=coli3, columnspan=2, sticky=(W, N, E, S))
    group3 = Label(mainframe4, text="Reaction Constant", padx=Pad_x, pady=Pad_y).grid(row=0, column=coli5, columnspan=2, sticky=(W, N, E, S))
    group4 = Label(mainframe4, text="Activation Energy", padx=Pad_x, pady=Pad_y).grid(row=0, column=coli7, columnspan=2, sticky=(W, N, E, S))
    Box_1 = Entry(mainframe4, width=7, textvariable=stringvarsr[i - 1])
    Box_1.grid(column=coli1, row=i, columnspan=CE, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    Label_0 = Label(mainframe4, text="Reaction {} ".format(i), padx=Pad_x, pady=Pad_y)
    Label_0.grid(column=coli0, row=i, sticky=(W, N, E, S))
    Label_0.rowconfigure(int(i), weight=1)
    Label_0.columnconfigure(coli0, weight=1)
    entriesr.append(Box_1)
    entriesp.append(Entry(mainframe4, width=7, textvariable=stringvarsp[i - 1]))
    entriesp[i - 1].grid(column=coli3, row=i, columnspan=CE, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    entriesk.append(Entry(mainframe4, width=7, textvariable=stringvarsk[i - 1]))
    entriesk[i - 1].grid(column=coli6, row=i, columnspan=1, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    if len(str(i)) >= 2:
        Label(mainframe4, text='k{}{}'.format(chr(0x2080 + int(str(i)[0])), chr(0x2080 + int(str(i)[-1])))).grid(column=coli5, row=i, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    elif len(str(i)) == 1:
        Label(mainframe4, text='k{}'.format(chr(0x2080 + int(str(i)[0])))).grid(column=coli5, row=i, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    entriesea.append(Entry(mainframe4, width=7, textvariable=stringvarsea[i - 1]))
    entriesea[i - 1].grid(column=coli8, row=i, columnspan=1, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    if len(str(i)) >= 2:
        Label(mainframe4, text='Ea{}{} [kJ/mol]'.format(chr(0x2080 + int(str(i)[0])), chr(0x2080 + int(str(i)[-1])))).grid(column=coli7, row=i, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    elif len(str(i)) == 1:
        Label(mainframe4, text='Ea{} [kJ/mol]'.format(chr(0x2080 + int(str(i)[0])))).grid(column=coli7, row=i, sticky=(W, N, E, S), padx=Pad_x, pady=Pad_y)
    entriesc.append(Checkbutton(mainframe4, text="Reaction {} Reversable".format(i), variable=intvars[i - 1]).grid(column=coli9, row=i, columnspan=2, sticky=(W, N, E, S)))
B4 = Button(root3, text="OK", command=close_window3).grid(column=coli9, row=1, padx=Pad_x, pady=Pad_y)
root3.mainloop()

root4 = Tk()
root4.title("Reactions")
mainframe4 = ttk.Frame(root4, padding="3 3 12 12")
mainframe4.grid(column=0, row=0, sticky=(N, W, E, S))
root4.columnconfigure(0, weight=1)
root4.rowconfigure(0, weight=1)
stringvarsr4 = []
stringvarsp4 = []
stringvarsrc = []
stringvarspc = []
entriesr4a = []
entriesp4a = []
entriesr4 = []
entriesp4 = []
entriesrca = []
entriespca = []
entriesrc = []
entriespc = []
rstrings = []

num_chems = len(chemnamesl)
for i in range(0, rxnnum, 1):
    rval = reverse[i]
    if rval == 0:
        rstrings.append(u"\u2192")
    elif rval == 1:
        rstrings.append(u"\u21CB")
for i in range(0, rxnnum, 1):
    stringvarsr4.append([StringVar(value="1") for i in range(0, reactants_num[i], 1)])
    stringvarsp4.append([StringVar(value="1") for i in range(0, products_num[i], 1)])
    stringvarsrc.append([StringVar(value="1") for i in range(0, reactants_num[i], 1)])
    stringvarspc.append([StringVar(value="1") for i in range(0, products_num[i], 1)])

for i in range(0, rxnnum, 1):
    mainframe4.rowconfigure(i + 1, weight=1)
    rangeval = int(reactants_num[i])
    rangep = int(products_num[i])
    rval2 = int(reactants_num[i])*2 + 1
    int1 = 1
    int2 = 2
    jval = 1
    for j in range(0, reactants_num[i], 1):
        mainframe4.columnconfigure(jval, weight=1)
        entriesr4a.append(Entry(mainframe4, width=7, textvariable=stringvarsr4[i][j]))
        entriesr4a[-1].grid(column=jval, row=int(i + 1))
        jval += 1
        mainframe4.columnconfigure(jval, weight=1)
        combbo = Combobox(mainframe4, values=chemnamesl)
        combbo.grid(column=jval, row=int(i + 1))
        jval += 1
        mainframe4.columnconfigure(jval, weight=1)
        entriesrca.append(combbo)
        if j < reactants_num[i]-1:
            mainframe4.columnconfigure(jval, weight=1)
            Label(mainframe4, text=" + ").grid(column=jval, row=int(i + 1))
            jval += 1
        elif j == reactants_num[i]-1:
            mainframe4.columnconfigure(jval, weight=1)
            Label(mainframe4, text=" {} ".format(rstrings[i])).grid(column=jval, row=int(i + 1))
            jval += 1
        int1 += 1
        int2 += 1
    endp = rval2 + products_num[i]+1
    entriesr4.append(entriesr4a[:])
    entriesr4a.clear()
    entriesrc.append(entriesrca[:])
    entriesrca.clear()
    for k in range(0, products_num[i], 1):
        int1b = 1
        int2b = 2
        mainframe4.columnconfigure(jval, weight=1)
        entriesp4a.append(Entry(mainframe4, width=7, textvariable=stringvarsp4[i][k]))
        entriesp4a[-1].grid(column=jval, row=int(i + 1))
        jval += 1
        mainframe4.columnconfigure(jval, weight=1)
        combbb = Combobox(mainframe4, values=chemnamesl)
        combbb.grid(column=jval, row=int(i + 1))
        jval += 1
        mainframe4.columnconfigure(jval, weight=1)
        entriespca.append(combbb)
        if k < products_num[i]-1:
            Label(mainframe4, text=" + ").grid(column=jval, row=int(i + 1))
            jval += 1
            mainframe4.columnconfigure(jval, weight=1)
        else:
            mainframe4.rowconfigure(int(rxnnum) + 2, weight=1)
            passB5 = Button(root4, text="OK", command=close_window4).grid(column=2, row=int(rxnnum+2))
    entriesp4.append(entriesp4a[:])
    entriesp4a.clear()
    entriespc.append(entriespca[:])
    entriespca.clear()
root4.mainloop()

rxns_strs = ["Reaction {}".format(int(i + 1)) for i in range(0, rxnnum, 1)]

for i in range(0, rxnnum, 1):
    indexnum = int(i + 1)
    keys = chemnamesl
    valuesr = coeffsr[i][:]
    valuesp = coeffsp[i][:]
    dictionary = {"Ea": indexnum, "K_Value": indexnum, "Reverse": reverse[i], "Reactants": dict(zip(keys, valuesr)), "Products": dict(zip(keys, valuesp))}
    Initreactions.append(dictionary)

for i in range(0, len(chemnamesl), 1):
    indexnum = int(i + 1)
    namev = chemnamesl[i]
    name_index = chemnamesl.index(namev)
    keys = rxns_strs
    valuesfor = [0*ij for ij in range(0, rxnnum, 1)]
    valuesrev = [0*ik for ik in range(0, rxnnum, 1)]
    for j in range(0, rxnnum, 1):
        valuef = coeffsr[j][name_index]
        if valuef != 0:
            valuesfor[j] = int(-1)
        elif coeffsp[j][name_index] != 0:
            valuesfor[j] = int(1)
        valuesrev[j] = coeffsp[j][name_index]*int(-1*reverse[j]) #coeffsp[j][name_index]*
    dictionary2 = {"Name": "{}".format(str(namev)), "Reactions": dict(zip(keys, valuesfor)), "Reverse": dict(zip(keys, valuesrev))}
    Eqlist.append(dictionary2)


class symbolgen:

    def __init__(self, nlist, Initlist, EQlist):
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
        return self.symfunc(self.nameslist, self.rxnnum)[0]

    def rinit(self):
        return self.initfunc(self.initreactions, self.C)[0]

    def pinit(self):
        return self.initfunc(self.initreactions, self.C)[1]

    initreactions = property(initl)
    C = property(symsinit, symsinit)
    reactants = property(rinit)
    products = property(pinit)
    latexs = property(latexin)

    def mfunci(funcsl, ylist, i, j):
        return diff(funcsl[i], ylist[j])

    def symfunc(names, rxnum):
        Csyms = [symbols(r'C_{}'.format('{}'.format(i))) for i in names]
        Ksyms = [symbols(r'K_{}'.format(j)) for j in range(rxnum)]
        EAsyms = [symbols(r'Ea_{}'.format(k)) for k in range(rxnum)]
        Tsyms = [symbols('T')]

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
                if A is True and count == 0:
                    As.append(vv)
                    SG.append(vv)
                    count = 0
                    fnum -= 1
                if A is True and count > 0:

                    NNa = "".join(Ns)
                    SG.append(NNa)
                    SG.append(vv)
                    Ns.clear()
                    count = 0
                    fnum -= 1
                if A is True and count >= 2:

                    NNa = "".join(Ns)
                    NNs.append(NNa)
                    Ns.clear()
                    SG.append(NNa)
                    SG.append(vv)
                    count = 0
                    fnum -= 1

                if N is True and ff > 1:
                    Ns.append(vv)
                    count += 1
                if N is True and ff <= 1:
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
                except Exception:
                    pass
            ss = "".join(SG)
            s.clear()
            fcs.append(ss)
        return fcs

    def rterm(Ci, a):
        termi = Mul(a, Pow(Ci, abs(int(a))))
        return termi

    def rprod(Ci, a, Cj, b):
        term1 = symbolgen.rterm(Ci, a)
        term2 = symbolgen.rterm(Cj, b)
        term3 = Mul(term1, term2)
        return term3

    def initfunc(initlist, C):
        reactants = []
        products = []
        for i, j in enumerate(initlist):
            Reactants = initlist[i]['Reactants']
            Products = initlist[i]['Products']
            Rvals = list(Reactants.values())
            Pvals = list(Products.values())
            Ks = symbols('k_{}'.format(i + 1))
            Eas = symbols('Ea_{}'.format(i + 1))
            RT = Mul(Symbol('R'), Symbol('T'))
            RTI = Pow(RT, Integer(-1))
            EART = Mul(Eas, RTI)
            EARTI = Mul(EART, Integer(-1))
            ee = exp(EARTI)
            rterms = []
            pterms = []
            rtotal = Integer(1)
            ptotal = Integer(1)
            for k, li in zip(C, Rvals):
                if li != 0:
                    term = symbolgen.rterm(k, li)
                    rterms.append(term)
            for t in rterms:
                rtotal = Mul(rtotal, t)
            for m, n in zip(C, Pvals):
                if n != 0:
                    pterm = symbolgen.rterm(m, n)
                    pterms.append(pterm)
            for tt in pterms:
                ptotal = Mul(ptotal, tt)
            reactants.append(Mul(Ks, Mul(rtotal, ee)))
            products.append(Mul(Ks, Mul(ptotal, ee)))
        return [reactants, products]

    def eqlist(eqlistl, R, P):
        reactants = R
        products = P
        EQS = []
        leqns = []
        for i, j in enumerate(eqlistl):
            Reactions = eqlistl[i]['Reactions']
            Reverse = eqlistl[i]['Reverse']
            Rxn = list(Reactions.values())
            RxnR = list(Reverse.values())
            eqn = []
            Reacts = [i*j for i, j in zip(Rxn, reactants) if i != 0]
            Prods = [i*j for i, j in zip(RxnR, products) if i != 0]
            if not Prods:
                eee = sum(Reacts)
                rlatex = latex(eee)
                leqns.append(rlatex)
                EQS.append(eee)
            else:
                eqn = sum(Reacts)
                peqn = sum(Prods)
                eeqn = Add(eqn, peqn)
                rlatex = latex(eeqn)
                leqns.append(rlatex)
                EQS.append(eeqn)
        return [EQS, leqns]

    def dislat(lnames, latexs, indvar):
        Latexs = []
        Displays = []
        for i in range(len(latexs)):
            dd = '{d' + 'C_{}'.format(symbols(lnames[i])) + '}'
            dt = '{d' + '{}'.format(symbols(indvar)) + '}'
            dde = r'${}\frac{}{}{}'.format(r"\mathbf{", dd, dt, "}") + ' = ' + '{}$'.format(latexs[i])
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
            rname = [symbols('{}'.format(Rname[h])) for h in range(len(Rname)) if Rvals[h] != 0]
            Pvals = list(Products.values())
            pvals = [Pvals[kk] for kk in range(len(Pvals)) if Pvals[kk] != 0]
            Pname = symbolgen.numfunc(list(Products.keys()))
            pname = [symbols('{}'.format(Pname[h])) for h in range(len(Pname)) if Pvals[h] != 0]
            CRvals = sum([Mul(Integer(ii), jj) for ii, jj in zip(rvals, rname) if ii != 0])
            CPvals = sum([Mul(Integer(ii), jj) for ii, jj in zip(pvals, pname) if ii != 0])
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
            eqn = [r'{}'.format(e).replace('{', '').replace('}', '')]
            Ksyms = [symbols('k_{}'.format(i + 1)) for i in range(len(kk))]
            EAsyms = [symbols('Ea_{}'.format(i + 1)) for i in range(len(ea))]
            kdictionary = dict(zip(Ksyms, kk))
            eadictionary = dict(zip(EAsyms, ea))
            eqn3 = e.subs(kdictionary)
            eqn4 = eqn3.subs(eadictionary)
            eqn5 = eqn4.subs({'R': 8.31446261815324})
            eqn6b = eqn5.subs({'*exp': '*sp.exp'})
            EQLISTF.append(eqn6b)
            EQLIST.append(eqn[0])

        return EQLIST, EQLISTF

    def jacobian(rhs, y):
        eqnl = len(rhs)
        cl = len(y)

        def mfunc(i, j):
            return diff(rhs[i], y[j])
        J = [[i for i in range(cl)] for j in range(eqnl)]
        Jf = [[sp.diff(rhs[j], y[i]) for i in range(cl)] for j in range(eqnl)]
        Jn = [[i for i in range(cl)] for j in range(eqnl)]
        Jm = [[i for i in range(cl)] for j in range(eqnl)]
        ix, jx = symbols("ix jx")
        Ja = Matrix(len(rhs), len(y), lambda i, j: mfunc(i, j))
        for i in range(eqnl):
            for j in range(cl):
                J[i][j] = str('{}'.format('{}'.format(mfunc(i, j)).replace('*exp', '*sp.exp')))
        for i in range(eqnl):
            for j in range(cl):
                Jn[i][j] = str('{}'.format('{}'.format(mfunc(i, j)).replace('*exp', '*np.exp')))
                Jm[i][j] = str('{}'.format('{}'.format(mfunc(i, j)).replace('*exp', '*math.exp')))
        MatrixJ = simplify(Matrix(Ja))
        LatexMatrix = sp.latex(matrix2numpy(Matrix(Jf)))
        lm = latex(MatrixJ, mode='inline', itex=True, mat_delim="(", mat_str='array')
        return J, Jn, Jm, MatrixJ, lm, LatexMatrix

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
        Cs.append("T")
        RHS, RHSf = symbolgen.rhseqs(equats, kk, ea, r)
        Jac, JacNumpy, JacMath, JacSimple, lm, latexmatrix = symbolgen.jacobian(RHSf, Cs)
        JacS, JacNumpyS, JacMathS, JacSimpleS, lmS, latexmatrixS = symbolgen.jacobian(RHS, Cs)
        symbolgen.csave(Chem, filepathf)
        symbolgen.psave(names, dlat, filepathf)
        KS = [str(r"{}".format(Ks[i])) for i in range(len(Ks))]
        EAS = [str(r"{}".format(EAs[i])) for i in range(len(EAs))]
        EAK = KS.copy()
        EAK.extend(EAS)
        symbolgen.fsave(filepathf, equats, dlat, Chem, ChemW, RHS, RHSf, Jac, JacNumpy, JacMath, JacSimple, lm, latexmatrix, JacS, JacNumpyS, JacMathS, JacSimpleS, lmS, latexmatrixS, Cs, EAK, names)
        return Cs, Ks, EAs, reacts, prods, equats, slat, dlat, Chem, ChemD, ChemW, RHS, RHSf, Jac, JacNumpy, JacMath, JacSimple, lm, latexmatrix, JacS, JacNumpyS, JacMathS, JacSimpleS, lmS, latexmatrixS

    def psave(nameslist, LATEXD, fpath):
        filename = fpath
        for s, k in enumerate(LATEXD):
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            left, width = .25, .5
            bottom, height = .25, .5
            right = left + width
            top = bottom + height
            ax.set_axis_off()
            ax.text(0.5 * (left + right), 0.5 * (bottom + top), k, va='center', ha='center', bbox=dict(boxstyle="round", fc="white", alpha=0.3, ec="black", pad=0.2))
            fig.savefig(r'{}\{}.svg'.format(filename, nameslist[s]), bbox_inches='tight')
            fig.savefig(r'{}\{}.pdf'.format(filename, nameslist[s]), bbox_inches='tight')
            fig.savefig(r'{}\{}.png'.format(filename, nameslist[s]), bbox_inches='tight')
            plt.close()

    def csave(LATEXC, fpath):
        filename = fpath
        for s, k in enumerate(LATEXC):
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
            left, width = .25, .5
            bottom, height = .25, .5
            right = left + width
            top = bottom + height
            ax.set_axis_off()
            ax.text(0.5*(left+right), 0.5*(bottom+top), k, va='center', ha='center', bbox=dict(boxstyle="round", fc="white", alpha=0.3, ec="black", pad=0.2))
            fig.savefig(r'{}\Reaction {}.svg'.format(filename, s + 1), bbox_inches='tight')
            fig.savefig(r'{}\Reaction {}.pdf'.format(filename, s + 1), bbox_inches='tight')
            fig.savefig(r'{}\Reaction {}.png'.format(filename, s + 1), bbox_inches='tight')
            plt.close()

    def fsave(ffpath, eqns, eqnslat, crxns, crxnsw, rhseq, rhseqf, Jac, JacN, JacMath, JacSimple, lm, latexmatrix, JacSy, JacSyN, JacMathSy, JacSimpleSy, lmSy, latexmatrixSy, C, EAK, nameslist):
        with open(r"{}\Equations.txt".format(ffpath), "w") as output:
            output.write("[")
            el = len(eqns)
            eel = 0
            for eqn in eqns:
                eel += 1
                if eel < el:
                    output.write('{},\n'.format(str(eqn)))
                if eel >= el:
                    output.write('{}]'.format(str(eqn)))

        with open(r"{}\EquationsLatex.txt".format(ffpath), "w") as output:
            for eqnlat in eqnslat:
                output.write('{}\n'.format(str(eqnlat)))

        with open(r"{}\EquationsTex.tex".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "$")
            output.write(r"\documentclass{article}")
            output.write("\n")
            output.write(r"\usepackage{amsmath, nccmath, bm}")
            output.write("\n")
            output.write(r"\usepackage[bottom=0.2in,top=0.2in,left=0.2in,right=0.2in]{geometry}")
            output.write("\n")
            output.write(r"\begin{document}")
            output.write("\n")
            output.write(r"\begin{fleqn}")
            output.write("\n")
            for eqnlat in eqnslat:
                output.write(r"\begin{equation}")
                output.write("\n")
                output.write(r"\begin{split}")
                output.write("\n")
                output.write('{}\n'.format(str(eqnlat).translate(removetable)))
                output.write(r"\end{split}")
                output.write("\n")
                output.write(r"\end{equation}")
                output.write("\n")
                output.write(r"\\")
                output.write("\n")
            output.write(r"\end{fleqn}")
            output.write("\n")
            output.write(r"\end{document}")

        with open(r"{}\ReactionsLatex.txt".format(ffpath), "w") as output:
            for crxn in crxns:
                output.write('{}\n'.format(str(crxn)))

        with open(r"{}\ReactionsLatexWord.txt".format(ffpath), "w") as output:
            for crxnw in crxnsw:
                output.write('{}\n'.format(str(crxnw)))

        with open(r"{}\RHSsymbols.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def RHS(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            ll = len(rhseq)
            eqsss = []
            for rhs in rhseq:
                lr = rhseq.index(rhs)
                if lr < ll:
                    eqsss.append(str("EQ_{}".format(nameslist[rhseq.index(rhs)])))
                    output.write("    EQ_{} = {}\n".format(nameslist[rhseq.index(rhs)], rhs))
                elif lr >= ll:
                    eqsss.append(str("EQ_{}".format(nameslist[rhseq.index(rhs)])))
                    output.write("    EQ_{} = {}\n".format(nameslist[rhseq.index(rhs)], rhs))
            output.write("    return [{}]".format(("{}".format(eqsss)).translate(removetableB)))

        with open(r"{}\RHS.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def RHS(t, y):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            ll = len(rhseqf)
            lr = 0
            eqsss = []
            for rhsff in rhseqf:
                lr += 1
                if lr < ll:
                    eqsss.append(str("EQ_{}".format(nameslist[rhseqf.index(rhsff)])))
                    output.write("    EQ_{} = {}\n".format(nameslist[rhseqf.index(rhsff)], rhsff))
                elif lr >= ll:
                    eqsss.append(str("EQ_{}".format(nameslist[rhseqf.index(rhsff)])))
                    output.write("    EQ_{} = {}\n".format(nameslist[rhseqf.index(rhsff)], rhsff))
            output.write("    return [{}]".format(("{}".format(eqsss)).translate(removetableB)))

        with open(r"{}\Jacobian.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacMathSy)
            jjj = 0
            for i in range(len(JacMathSy)):
                jjj += 1
                Jrow = JacMathSy[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianSympy.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacSy)
            jjj = 0
            for i in range(len(JacSy)):
                jjj += 1
                Jrow = JacSy[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianNumpy.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacN)
            jjj = 0
            for i in range(len(JacN)):
                jjj += 1
                Jrow = JacN[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianMatrix.txt".format(ffpath), 'w') as output:
            output.write('{}'.format(JacSimple))

        with open(r"{}\JacobianLatex.txt".format(ffpath), "w") as output:
            output.write('{}'.format(lm))

        with open(r"{}\RHS.txt".format(ffpath)) as filein, open(r"{}\RightHandSide.txt".format(ffpath), 'w') as fileout:
            fileinl = filein.readlines()
            lfia = len(fileinl)
            lffb = 0
            for line in fileinl:
                lffb += 1
                line = line.replace("'", "")
                line = line.replace("exp", "sp.exp")
                if lffb < lfia:
                    fileout.write('{}'.format(line))
                elif lffb >= lfia:
                    fileout.write('{}'.format(line))

        with open(r"{}\RHSsymbols.txt".format(ffpath)) as filein, open(r"{}\RightHandSideSymbols.txt".format(ffpath), 'w') as fileout:
            fileinl = filein.readlines()
            lfi = len(fileinl)
            lff = 0
            for line in fileinl:
                line = line.replace("'", "")
                line = line.replace("exp", "math.exp")
                lff += 1
                if lff < lfi:
                    fileout.write('{}'.format(line))
                elif lff >= lfi:
                    fileout.write('{}'.format(line))

        pickle.dumps(JacSimple)
        with open(r'{}\JacobianMatrixPickle.txt'.format(ffpath), 'wb') as f:
            pickle.dump(JacSimple, f)

        with open(r"{}\JacobianSymbolic.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacMathSy)
            jjj = 0
            for i in range(len(JacMathSy)):
                jjj += 1
                Jrow = JacMathSy[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianSymbolicSympy.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacSy)
            jjj = 0
            for i in range(len(JacSy)):
                jjj += 1
                Jrow = JacSy[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianSymbolicNumpy.txt".format(ffpath), "w") as output:
            removetable = str.maketrans('', '', "[]'")
            removetableB = str.maketrans('', '', "[]'")
            output.write("def Jacob(t, y, *args):\n")
            output.write("    {} = args\n".format(str("{}".format(EAK)).translate(removetable)))
            output.write("    {} = y\n".format(str("{}".format(C)).translate(removetable)))
            output.write("    Jac = [")
            jj = len(JacSy)
            jjj = 0
            for i in range(len(JacSyN)):
                jjj += 1
                Jrow = JacSyN[i][:]
                if i == 0:
                    output.write(('{},\n'.format(Jrow)).replace("'", ""))
                if jjj < jj and i != 0:
                    output.write(('           {},\n'.format(Jrow)).replace("'", ""))
                elif jjj >= jj:
                    output.write(('           {}'.format(Jrow)).replace("'", ""))
            output.write("]\n")
            output.write("    return Jac")

        with open(r"{}\JacobianMatrixSymbolic.txt".format(ffpath), 'w') as output:
            output.write('{}'.format(JacSimpleSy))

        with open(r"{}\JacobianLatexSymbolic.txt".format(ffpath), "w") as output:
            output.write('{}'.format(lmSy))

        pickle.dumps(JacSimpleSy)
        with open(r'{}\JacobianMatrixPickleSymbolic.txt'.format(ffpath), 'wb') as f:
            pickle.dump(JacSimpleSy, f)

        try:
            create_pdf(r"{}\EquationsTex.tex".format(ffpath), "{}".format(ffpath))
        except Exception:
            print("Coulnd't convert tex file.")
            pass


ea = [i*1000.0 for i in eaf]
RRv = RR[0]
rxnnumf = rxnsvl[0]
C, KKS, EAS, reacts, prods, equations, slat, dlat, chem, chemD, chemw, rhs, rhsf, jac, jacnumpy, Jacmath, JacSimple, lm, latexmatrix, jacsy, jacnumpysy, jacmathsy, jacsimplesy, lmsy, latexmatrixsy = symbolgen.fullgen(chemnamesl, rxnnumf, Initreactions, Eqlist, indvdf[0], ffpath[0], kk, ea, RRv)
# print(rhs)
# print(rhsf)
