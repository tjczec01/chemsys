# Symbol Generator

*symbolgen.py* is a Gui based chemical reaction system generator that generates the right hand side (**RHS**) of a chemical system to be used with the [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint)  method, the **jacobian**, and the overall equations for each reations and each individual chemical species. It will generate **Latex** formatted equations for the individual chemical reactions, overall mass balances for each species, and both symbolic and numerical text files for the **RHS** and **jacobian**.


# 1. Chemical reactions

This program will generate the **Latex** forms of each individual reaction and then save them as both a **pdf** and **svg**. The string forms of the equations will be saved in a **text** (txt) file. An example of some inital reactions are given below.

## <center> <b>List of Reactions</b>
    
$$EDC \longrightarrow R_{1} + R_{2}$$
$$EDC + R_{1} \longrightarrow HCl + R_{3}$$
$$EDC + R_{2} \longrightarrow EC + R_{3}$$
$$EDC + R_{4} \longrightarrow C_{11}$$
$$EDC + R_{5} \longrightarrow R_{3} + VCM$$
$$EDC + R_{6} \longrightarrow C_{112} + R_{3}$$
$$R_{1} + R_{2} \longrightarrow HCl + VCM$$
$$R_{1} + R_{3} \longrightarrow Di + HCl$$
$$EC + R_{1} \longrightarrow HCl + R_{2}$$
$$C_{11} + R_{1} \longrightarrow HCl + R_{4}$$
$$C_{112} + R_{1} \longrightarrow HCl + R_{6}$$
$$R_{1} + VCM \longrightarrow R_{4}$$
$$R_{1} + VCM \longrightarrow HCl + R_{5}$$
$$R_{2} + VCM \longrightarrow EC + R_{5}$$
$$R_{4} + VCM \longrightarrow C_{4}H_{6}Cl_{2} + R_{1}$$
$$R_{5} + VCM \longrightarrow CP + R_{1}$$
$$R_{3} \rightleftharpoons R_{1} + VCM$$
$$R_{5} \rightleftharpoons C_{2}H_{2} + R_{1}$$
$$R_{6} \rightleftharpoons Di + R_{1}$$
$$2*C_{2}H_{2} + R_{5} \longrightarrow C_{6}H_{6} + R_{1}$$
$C_{2}H_{2} + 2*R_{1} \longrightarrow 2*Coke + 2*HCl$

# 2. Overall reactions for each chemical species

This program will generate the **Latex** forms of each individual reaction and then save them as both a **pdf** and **svg**. The string forms of the equations will be saved in a **text** (txt) file. An example of some inital reactions are given below.



<left> $\dfrac{\textbf{dEDC}}{\textbf{dz}} = - C_{EDC} C_{R1} K_{2} e^{\frac{Ea_{2}}{R T}} - C_{EDC} C_{R2} K_{3} e^{\frac{Ea_{3}}{R T}} - C_{EDC} C_{R4} K_{4} e^{\frac{Ea_{4}}{R T}} - C_{EDC} C_{R5} K_{5} e^{\frac{Ea_{5}}{R T}} - C_{EDC} C_{R6} K_{6} e^{\frac{Ea_{6}}{R T}} - C_{EDC} K_{1} e^{\frac{Ea_{1}}{R T}}$






# 3. Right Hand Side 


The right hand side (**RHS**) of the system of equations will be generated both symbolically and with initial values substituted into their respective places. This is the required system for scipy's [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint) method. This function requires a callable in the form of 

```python
def RHS(t, y):

    return [-C_EDC*C_R1*K_2*sp.exp(Ea_2/(R*T)) - C_EDC*C_R2*K_3*sp.exp(Ea_3/(R*T)) - C_EDC*C_R4*K_4*sp.exp(Ea_4/(R*T)) - C_EDC*C_R5*K_5*sp.exp(Ea_5/(R*T)) - C_EDC*C_R6*K_6*sp.exp(Ea_6/(R*T)) - C_EDC*K_1*sp.exp(Ea_1/(R*T)),
            -C_EC*C_R1*K_9*sp.exp(Ea_9/(R*T)) + C_EDC*C_R2*K_3*sp.exp(Ea_3/(R*T)) + C_R2*C_VCM*K_14*sp.exp(Ea_14/(R*T)),
            C_C11*C_R1*K_10*sp.exp(Ea_10/(R*T)) + C_C112*C_R1*K_11*sp.exp(Ea_11/(R*T)) + 2*C_C2H2*C_R1**2*K_21*sp.exp(Ea_21/(R*T)) + C_EC*C_R1*K_9*sp.exp(Ea_9/(R*T)) + C_EDC*C_R1*K_2*sp.exp(Ea_2/(R*T)) + C_R1*C_R2*K_7*sp.exp(Ea_7/(R*T)) + C_R1*C_R3*K_8*sp.exp(Ea_8/(R*T)) + C_R1*C_VCM*K_13*sp.exp(Ea_13/(R*T)),
            2*C_C2H2*C_R1**2*K_21*sp.exp(Ea_21/(R*T)),
            C_R5*C_VCM*K_16*sp.exp(Ea_16/(R*T)),
            -C_Di*C_R1*K_19*sp.exp(Ea_19/(R*T)) + C_R1*C_R3*K_8*sp.exp(Ea_8/(R*T)) + C_R6*K_19*sp.exp(Ea_19/(R*T)),
            C_R4*C_VCM*K_15*sp.exp(Ea_15/(R*T)),
            2*C_C2H2**2*C_R5*K_20*sp.exp(Ea_20/(R*T)),
            -2*C_C2H2**2*C_R5*K_20*sp.exp(Ea_20/(R*T)) - 2*C_C2H2*C_R1**2*K_21*sp.exp(Ea_21/(R*T)) - C_C2H2*C_R1*K_18*sp.exp(Ea_18/(R*T)) - C_R5*K_18*sp.exp(Ea_18/(R*T)),
            -C_C11*C_R1*K_10*sp.exp(Ea_10/(R*T)) + C_EDC*C_R4*K_4*sp.exp(Ea_4/(R*T)),
            -C_C112*C_R1*K_11*sp.exp(Ea_11/(R*T)) + C_EDC*C_R5*K_5*sp.exp(Ea_5/(R*T)),
            -C_C11*C_R1*K_10*sp.exp(Ea_10/(R*T)) - C_C112*C_R1*K_11*sp.exp(Ea_11/(R*T)) + 2*C_C2H2**2*C_R5*K_20*sp.exp(Ea_20/(R*T)) - 2*C_C2H2*C_R1**2*K_21*sp.exp(Ea_21/(R*T)) - C_C2H2*C_R1*K_18*sp.exp(Ea_18/(R*T)) - C_Di*C_R1*K_19*sp.exp(Ea_19/(R*T)) - C_EC*C_R1*K_9*sp.exp(Ea_9/(R*T)) - C_EDC*C_R1*K_2*sp.exp(Ea_2/(R*T)) + C_EDC*K_1*sp.exp(Ea_1/(R*T)) - C_R1*C_R2*K_7*sp.exp(Ea_7/(R*T)) - C_R1*C_R3*K_8*sp.exp(Ea_8/(R*T)) - C_R1*C_VCM*K_12*sp.exp(Ea_12/(R*T)) - C_R1*C_VCM*K_13*sp.exp(Ea_13/(R*T)) - C_R1*C_VCM*K_17*sp.exp(Ea_17/(R*T)) + C_R3*K_17*sp.exp(Ea_17/(R*T)) + C_R4*C_VCM*K_15*sp.exp(Ea_15/(R*T)) + C_R5*C_VCM*K_16*sp.exp(Ea_16/(R*T)) + C_R5*K_18*sp.exp(Ea_18/(R*T)) + C_R6*K_19*sp.exp(Ea_19/(R*T)),
            C_EC*C_R1*K_9*sp.exp(Ea_9/(R*T)) - C_EDC*C_R2*K_3*sp.exp(Ea_3/(R*T)) + C_EDC*K_1*sp.exp(Ea_1/(R*T)) - C_R1*C_R2*K_7*sp.exp(Ea_7/(R*T)) - C_R2*C_VCM*K_14*sp.exp(Ea_14/(R*T)),
            C_EDC*C_R1*K_2*sp.exp(Ea_2/(R*T)) + C_EDC*C_R2*K_3*sp.exp(Ea_3/(R*T)) + C_EDC*C_R4*K_4*sp.exp(Ea_4/(R*T)) + C_EDC*C_R5*K_5*sp.exp(Ea_5/(R*T)) + C_EDC*C_R6*K_6*sp.exp(Ea_6/(R*T)) - C_R1*C_R3*K_8*sp.exp(Ea_8/(R*T)) + C_R1*C_VCM*K_17*sp.exp(Ea_17/(R*T)) - C_R3*K_17*sp.exp(Ea_17/(R*T)),
            C_C11*C_R1*K_10*sp.exp(Ea_10/(R*T)) - C_EDC*C_R4*K_4*sp.exp(Ea_4/(R*T)) + C_R1*C_VCM*K_12*sp.exp(Ea_12/(R*T)) - C_R4*C_VCM*K_15*sp.exp(Ea_15/(R*T)),
            2*C_C2H2**2*C_R5*K_20*sp.exp(Ea_20/(R*T)) + C_C2H2*C_R1*K_18*sp.exp(Ea_18/(R*T)) + C_EDC*C_R5*K_5*sp.exp(Ea_5/(R*T)) - C_R1*C_VCM*K_13*sp.exp(Ea_13/(R*T)) - C_R2*C_VCM*K_14*sp.exp(Ea_14/(R*T)) + C_R5*C_VCM*K_16*sp.exp(Ea_16/(R*T)) - C_R5*K_18*sp.exp(Ea_18/(R*T)),
            C_C112*C_R1*K_11*sp.exp(Ea_11/(R*T)) + C_Di*C_R1*K_19*sp.exp(Ea_19/(R*T)) - C_EDC*C_R6*K_6*sp.exp(Ea_6/(R*T)) - C_R6*K_19*sp.exp(Ea_19/(R*T)),
            C_EDC*C_R5*K_5*sp.exp(Ea_5/(R*T)) + C_R1*C_R2*K_7*sp.exp(Ea_7/(R*T)) - C_R1*C_VCM*K_12*sp.exp(Ea_12/(R*T)) - C_R1*C_VCM*K_13*sp.exp(Ea_13/(R*T)) - C_R1*C_VCM*K_17*sp.exp(Ea_17/(R*T)) - C_R2*C_VCM*K_14*sp.exp(Ea_14/(R*T)) + C_R3*K_17*sp.exp(Ea_17/(R*T)) - C_R4*C_VCM*K_15*sp.exp(Ea_15/(R*T)) - C_R5*C_VCM*K_16*sp.exp(Ea_16/(R*T))]
```


# 4. Jacobian 

The **Jacobian matrix** is symbolically generated in order to improve the accuracy and speed of the solvers used in the aforementioned method.



$\mathbf{f(t, y)} =\mathbf{f(t, C_{EDC}, C_{EC}, C_{HCl}, C_{Coke}, C_{CP}, C_{Di}, C_{C4H6Cl2}, C_{C6H6}, C_{C2H2}, C_{C11}, C_{C112}, C_{R1}, C_{R2}, C_{R3}, C_{R4}, C_{R5}, C_{R6}, C_{VCM}, T)}$






$$\mathbf{J}
=
\frac{d \mathbf{f}}{d \mathbf{y}}
=
\left[ \frac{\partial \mathbf{f}}{\partial y_1}
\cdots \frac{\partial \mathbf{f}}{\partial y_n} \right] 
=
\begin{bmatrix}
\frac{\partial f_1}{\partial y_1} & \cdots &
\frac{\partial f_1}{\partial y_n} \\
\vdots & \ddots & \vdots \\
\frac{\partial f_m}{\partial y_1} & \cdots & 
\frac{\partial f_m}{\partial y_n}
\end{bmatrix}$$


```python
def jacob(t, y):
    
     JJ = [[-13000000000000.0*C_R1*sp.exp(841.906485299082/T) - 1000000000000.0*C_R2*sp.exp(5051.43891179449/T) - 500000000000.0*C_R4*sp.exp(5412.25597692267/T) - 12000000000000.0*C_R5*sp.exp(4089.26007145269/T) - 200000000000.0*C_R6*sp.exp(5773.07304205085/T) - 5.9e+15*sp.exp(41133.1454246123/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -13000000000000.0*C_EDC*sp.exp(841.906485299082/T), -1000000000000.0*C_EDC*sp.exp(5051.43891179449/T), 0, -500000000000.0*C_EDC*sp.exp(5412.25597692267/T), -12000000000000.0*C_EDC*sp.exp(4089.26007145269/T), -200000000000.0*C_EDC*sp.exp(5773.07304205085/T), 0],
            [1000000000000.0*C_R2*sp.exp(5051.43891179449/T), -17000000000000.0*C_R1*sp.exp(481.089420170904/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, -17000000000000.0*C_EC*sp.exp(481.089420170904/T), 1000000000000.0*C_EDC*sp.exp(5051.43891179449/T) + 500000000000.0*C_VCM*sp.exp(3728.44300632451/T), 0, 0, 0, 0, 500000000000.0*C_R2*sp.exp(3728.44300632451/T)],
            [13000000000000.0*C_R1*sp.exp(841.906485299082/T), 17000000000000.0*C_R1*sp.exp(481.089420170904/T), 0, 0, 0, 0, 0, 0, 320000000000000.0*C_R1**2*sp.exp(8419.06485299082/T), 12000000000000.0*C_R1*sp.exp(721.634130256356/T), 17000000000000.0*C_R1*sp.exp(1804.08532564089/T), 12000000000000.0*C_C11*sp.exp(721.634130256356/T) + 17000000000000.0*C_C112*sp.exp(1804.08532564089/T) + 640000000000000.0*C_C2H2*C_R1*sp.exp(8419.06485299082/T) + 17000000000000.0*C_EC*sp.exp(481.089420170904/T) + 13000000000000.0*C_EDC*sp.exp(841.906485299082/T) + 10000000000000.0*C_R2*sp.exp(1563.54061555544/T) + 10000000000000.0*C_R3*sp.exp(1443.26826051271/T) + 120000000000000.0*C_VCM*sp.exp(6735.25188239266/T), 10000000000000.0*C_R1*sp.exp(1563.54061555544/T), 10000000000000.0*C_R1*sp.exp(1443.26826051271/T), 0, 0, 0, 120000000000000.0*C_R1*sp.exp(6735.25188239266/T)],
            [0, 0, 0, 0, 0, 0, 0, 0, 320000000000000.0*C_R1**2*sp.exp(8419.06485299082/T), 0, 0, 640000000000000.0*C_C2H2*C_R1*sp.exp(8419.06485299082/T), 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 300000000000.0*C_VCM*sp.exp(7336.61365760629/T), 0, 300000000000.0*C_R5*sp.exp(7336.61365760629/T)],
            [0, 0, 0, 0, 0, -20000000000000.0*C_R1*sp.exp(8419.06485299082/T), 0, 0, 0, 0, 0, -20000000000000.0*C_Di*sp.exp(8419.06485299082/T) + 10000000000000.0*C_R3*sp.exp(1443.26826051271/T), 0, 10000000000000.0*C_R1*sp.exp(1443.26826051271/T), 0, 0, 20000000000000.0*sp.exp(8419.06485299082/T), 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20000000000.0*C_VCM*sp.exp(3608.17065128178/T), 0, 0, 20000000000.0*C_R4*sp.exp(3608.17065128178/T)],
            [0, 0, 0, 0, 0, 0, 0, 0, 400000000000000.0*C_C2H2*C_R5*sp.exp(2405.44710085452/T), 0, 0, 0, 0, 0, 0, 200000000000000.0*C_C2H2**2*sp.exp(2405.44710085452/T), 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -400000000000000.0*C_C2H2*C_R5*sp.exp(2405.44710085452/T) - 320000000000000.0*C_R1**2*sp.exp(8419.06485299082/T) - 500000000000000.0*C_R1*sp.exp(10824.5119538453/T), 0, 0, -640000000000000.0*C_C2H2*C_R1*sp.exp(8419.06485299082/T) - 500000000000000.0*C_C2H2*sp.exp(10824.5119538453/T), 0, 0, 0, -200000000000000.0*C_C2H2**2*sp.exp(2405.44710085452/T) - 500000000000000.0*sp.exp(10824.5119538453/T), 0, 0],
            [500000000000.0*C_R4*sp.exp(5412.25597692267/T), 0, 0, 0, 0, 0, 0, 0, 0, -12000000000000.0*C_R1*sp.exp(721.634130256356/T), 0, -12000000000000.0*C_C11*sp.exp(721.634130256356/T), 0, 0, 500000000000.0*C_EDC*sp.exp(5412.25597692267/T), 0, 0, 0],
            [12000000000000.0*C_R5*sp.exp(4089.26007145269/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, -17000000000000.0*C_R1*sp.exp(1804.08532564089/T), -17000000000000.0*C_C112*sp.exp(1804.08532564089/T), 0, 0, 0, 12000000000000.0*C_EDC*sp.exp(4089.26007145269/T), 0, 0],
            [-13000000000000.0*C_R1*sp.exp(841.906485299082/T) + 5.9e+15*sp.exp(41133.1454246123/T), -17000000000000.0*C_R1*sp.exp(481.089420170904/T), 0, 0, 0, -20000000000000.0*C_R1*sp.exp(8419.06485299082/T), 0, 0, 400000000000000.0*C_C2H2*C_R5*sp.exp(2405.44710085452/T) - 320000000000000.0*C_R1**2*sp.exp(8419.06485299082/T) - 500000000000000.0*C_R1*sp.exp(10824.5119538453/T), -12000000000000.0*C_R1*sp.exp(721.634130256356/T), -17000000000000.0*C_R1*sp.exp(1804.08532564089/T), -12000000000000.0*C_C11*sp.exp(721.634130256356/T) - 17000000000000.0*C_C112*sp.exp(1804.08532564089/T) - 640000000000000.0*C_C2H2*C_R1*sp.exp(8419.06485299082/T) - 500000000000000.0*C_C2H2*sp.exp(10824.5119538453/T) - 20000000000000.0*C_Di*sp.exp(8419.06485299082/T) - 17000000000000.0*C_EC*sp.exp(481.089420170904/T) - 13000000000000.0*C_EDC*sp.exp(841.906485299082/T) - 10000000000000.0*C_R2*sp.exp(1563.54061555544/T) - 10000000000000.0*C_R3*sp.exp(1443.26826051271/T) - 120000000000000.0*C_VCM*sp.exp(6735.25188239266/T) - 210000000000000.0*C_VCM*sp.exp(10102.877823589/T) - 91000000000.0*C_VCM, -10000000000000.0*C_R1*sp.exp(1563.54061555544/T), -10000000000000.0*C_R1*sp.exp(1443.26826051271/T) + 210000000000000.0*sp.exp(10102.877823589/T), 20000000000.0*C_VCM*sp.exp(3608.17065128178/T), 200000000000000.0*C_C2H2**2*sp.exp(2405.44710085452/T) + 300000000000.0*C_VCM*sp.exp(7336.61365760629/T) + 500000000000000.0*sp.exp(10824.5119538453/T), 20000000000000.0*sp.exp(8419.06485299082/T), -120000000000000.0*C_R1*sp.exp(6735.25188239266/T) - 210000000000000.0*C_R1*sp.exp(10102.877823589/T) - 91000000000.0*C_R1 + 20000000000.0*C_R4*sp.exp(3608.17065128178/T) + 300000000000.0*C_R5*sp.exp(7336.61365760629/T)],
            [-1000000000000.0*C_R2*sp.exp(5051.43891179449/T) + 5.9e+15*sp.exp(41133.1454246123/T), 17000000000000.0*C_R1*sp.exp(481.089420170904/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 17000000000000.0*C_EC*sp.exp(481.089420170904/T) - 10000000000000.0*C_R2*sp.exp(1563.54061555544/T), -1000000000000.0*C_EDC*sp.exp(5051.43891179449/T) - 10000000000000.0*C_R1*sp.exp(1563.54061555544/T) - 500000000000.0*C_VCM*sp.exp(3728.44300632451/T), 0, 0, 0, 0, -500000000000.0*C_R2*sp.exp(3728.44300632451/T)],
            [13000000000000.0*C_R1*sp.exp(841.906485299082/T) + 1000000000000.0*C_R2*sp.exp(5051.43891179449/T) + 500000000000.0*C_R4*sp.exp(5412.25597692267/T) + 12000000000000.0*C_R5*sp.exp(4089.26007145269/T) + 200000000000.0*C_R6*sp.exp(5773.07304205085/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13000000000000.0*C_EDC*sp.exp(841.906485299082/T) - 10000000000000.0*C_R3*sp.exp(1443.26826051271/T) + 210000000000000.0*C_VCM*sp.exp(10102.877823589/T), 1000000000000.0*C_EDC*sp.exp(5051.43891179449/T), -10000000000000.0*C_R1*sp.exp(1443.26826051271/T) - 210000000000000.0*sp.exp(10102.877823589/T), 500000000000.0*C_EDC*sp.exp(5412.25597692267/T), 12000000000000.0*C_EDC*sp.exp(4089.26007145269/T), 200000000000.0*C_EDC*sp.exp(5773.07304205085/T), 210000000000000.0*C_R1*sp.exp(10102.877823589/T)],
            [-500000000000.0*C_R4*sp.exp(5412.25597692267/T), 0, 0, 0, 0, 0, 0, 0, 0, 12000000000000.0*C_R1*sp.exp(721.634130256356/T), 0, 12000000000000.0*C_C11*sp.exp(721.634130256356/T) + 91000000000.0*C_VCM, 0, 0, -500000000000.0*C_EDC*sp.exp(5412.25597692267/T) - 20000000000.0*C_VCM*sp.exp(3608.17065128178/T), 0, 0, 91000000000.0*C_R1 - 20000000000.0*C_R4*sp.exp(3608.17065128178/T)],
            [12000000000000.0*C_R5*sp.exp(4089.26007145269/T), 0, 0, 0, 0, 0, 0, 0, 400000000000000.0*C_C2H2*C_R5*sp.exp(2405.44710085452/T) + 500000000000000.0*C_R1*sp.exp(10824.5119538453/T), 0, 0, 500000000000000.0*C_C2H2*sp.exp(10824.5119538453/T) - 120000000000000.0*C_VCM*sp.exp(6735.25188239266/T), -500000000000.0*C_VCM*sp.exp(3728.44300632451/T), 0, 0, 200000000000000.0*C_C2H2**2*sp.exp(2405.44710085452/T) + 12000000000000.0*C_EDC*sp.exp(4089.26007145269/T) + 300000000000.0*C_VCM*sp.exp(7336.61365760629/T) - 500000000000000.0*sp.exp(10824.5119538453/T), 0, -120000000000000.0*C_R1*sp.exp(6735.25188239266/T) - 500000000000.0*C_R2*sp.exp(3728.44300632451/T) + 300000000000.0*C_R5*sp.exp(7336.61365760629/T)],
            [-200000000000.0*C_R6*sp.exp(5773.07304205085/T), 0, 0, 0, 0, 20000000000000.0*C_R1*sp.exp(8419.06485299082/T), 0, 0, 0, 0, 17000000000000.0*C_R1*sp.exp(1804.08532564089/T), 17000000000000.0*C_C112*sp.exp(1804.08532564089/T) + 20000000000000.0*C_Di*sp.exp(8419.06485299082/T), 0, 0, 0, 0, -200000000000.0*C_EDC*sp.exp(5773.07304205085/T) - 20000000000000.0*sp.exp(8419.06485299082/T), 0][12000000000000.0*C_R5*sp.exp(4089.26007145269/T), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10000000000000.0*C_R2*sp.exp(1563.54061555544/T) - 120000000000000.0*C_VCM*sp.exp(6735.25188239266/T) - 210000000000000.0*C_VCM*sp.exp(10102.877823589/T) - 91000000000.0*C_VCM, 10000000000000.0*C_R1*sp.exp(1563.54061555544/T) - 500000000000.0*C_VCM*sp.exp(3728.44300632451/T), 210000000000000.0*sp.exp(10102.877823589/T), -20000000000.0*C_VCM*sp.exp(3608.17065128178/T), 12000000000000.0*C_EDC*sp.exp(4089.26007145269/T) - 300000000000.0*C_VCM*sp.exp(7336.61365760629/T), 0, -120000000000000.0*C_R1*sp.exp(6735.25188239266/T) - 210000000000000.0*C_R1*sp.exp(10102.877823589/T) - 91000000000.0*C_R1 - 500000000000.0*C_R2*sp.exp(3728.44300632451/T) - 20000000000.0*C_R4*sp.exp(3608.17065128178/T) - 300000000000.0*C_R5*sp.exp(7336.61365760629/T)]]


     return sp.lambdify((t, y), JJ)
```


```python

```
