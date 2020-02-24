
# Class symbolgen
## List and description of functions

### 1. symbolgen.symfunc(nameslist, number_of_reactions)
  #### Creates a list of symbols given names as a list of strings.
  ```nameslist = ['EDC','EC','HCl','Coke', 'CP','Di','C4H6Cl2','C6H6','C2H2','C11','C112','R1','R2','R3','R4','R5','R6','VCM']
     number_of_reactions = len(Initreactions)
     Csymbols, K, EA, T = symfunc(nameslist, number_of_reactions)
  ```
### 2. symbolgen.init(Initreactions, Csymbols)
  #### Generates two list. The first is a list of all indivisual reactants for each reaction. The second is a list of all products for each reaction. 
  ```
  reactants, products = init(Initreactions, Csymbols)
  ``` 

### 3. symbolgen.eqlist(Eqlist, reactants, products)
 #### Uses the products and reactants generated above to product a list of overall equations for each individual chemical species. 
 ```
 equations, latexs = eqlist(Eqlist, reactants, products)
 ``` 
### 4. symbolgen.dislat(nameslist, latexs, indvar)
#### Generates the Latex/latex forms of the equations that are used to display and save the equations to pdf/svg. 
```
slatex, dlatex = dislat(nameslist, latexs, indvar)
``` 

### 5. symbolgen.rhseqs(equations)
#### Transforms the equations produced above into a more usefull form to be used with the solve_ivp method
```
reqs = rhseqs(equations)
``` 


