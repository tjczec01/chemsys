# Chemical System Generator
Generator for symbolic functions for both full and simple systems
* [Symbolic Generator](https://github.com/tjczec01/chemsys/blob/master/symbolgen.ipynb)

# Symbol Generator

**chemsys.py** is a Gui based chemical reaction system generator that generates the right hand side (**RHS**) of a chemical system to be used with the [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint)  method, the **jacobian**, and the overall equations for each reactions and each individual chemical species. It will generate **Latex** formatted equations for the individual chemical reactions, overall mass balances for each species, and both symbolic and numerical text files for the **RHS** and **jacobian**.

# Installation 
  ```pip install chemsys```  
  
 ```conda install -c tjczec01 chemsys``` 


# 1. Chemical reactions

This program will generate the **Latex** forms of each individual reaction and then save them as both a **pdf** and **svg**. The string forms of the equations will be saved in a **text** (txt) file. An example of some inital reactions are given below.
* [Symbolic Generator](https://github.com/tjczec01/symbolgen/blob/master/symbolgen.ipynb)

# 2. Overall reactions for each chemical species

This program will generate the **Latex** forms of each individual reaction and then save them as both a **pdf** and **svg**. The string forms of the equations will be saved in a **text** (txt) file. An example of some inital reactions are given below.
* [Symbolic Generator](https://github.com/tjczec01/symbolgen/blob/master/symbolgen.ipynb)

# 3. Right Hand Side 


The right hand side (**RHS**) of the system of equations will be generated both symbolically and with initial values substituted into their respective places. This is the required system for scipy's [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint) method. This function requires a callable in the form of  ``` fun(t, y):```
* [Symbolic Generator](https://github.com/tjczec01/symbolgen/blob/master/symbolgen.ipynb)

# 4. Jacobian 

The **Jacobian matrix** is symbolically generated in order to improve the accuracy and speed of the solvers used in the aforementioned method.
* [Symbolic Generator](https://github.com/tjczec01/symbolgen/blob/master/symbolgen.ipynb)

# 5. Use

# 5.1 Generate all necessary lists and values as follows.

```
from chemsys.cs import gui, symbolgen
chemical_names, number_of_reactions, Initial_reactions, Equation_list, indvdf, filepath, kvalues, ea_values, r_gas, pdf_path = gui.fullgui()
```

# 5.2 Calculate the jacobian and all other desired functions as follows.

```
C_Symbols, K_Vals, EA_Vals, reactants, products, equations, slat, dlat, chem, chemD, chemw, right_hand_side, right_hand_sidef, Jac, Jac_numpy, Jac_math, Jac_Simple, lm, latex_matrix, Jac_symbol, Jac_numpy_symbol, Jac_math_symbol, Jac_simple_symbol, lm_symbol, latex_matrix_symbol, DLatb = symbolgen.fullgen(chemical_names, number_of_reactions, Initial_reactions, Equation_list, indvdf, filepath, kvalues, ea_values, r_gas, chemical_names, pdf_path)

```

# 6.1 References

Czechorski, Travis J, "A kinetic and thermodynamic model of ethylene dichloride pyrolysis." (2019). Electronic Theses and Dissertations. Paper 3359.
Retrieved from https://ir.library.louisville.edu/etd/3359 

Reference the thesis by Travis Czechorski under the supervision of Dr. Vance Jaeger.
