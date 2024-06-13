
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Lagrange Methode")

st.write('The following page shows the code for a N-fold Pendulum calculated with the Lagrange Methode. This program is the generalization for the Double Pendulum. In the next section I will show the code and afterwards the results. The code itself is commented but the functions themself left out for better visibility. Documentation for the functions can be found under "Documentation Lagrange Methode". I will show in an alternating sequence code snippets (as implemented in jupyter notebook) and the corresponding output. The code itself can be found in the github repository linked in the "Homepage".')

st.write('Before the code and the results themself some remarks:')

st.write('1. I excluded the code for the functions in the snippets below for better visibility. The code for the functions is still inserted in the jupyter notebooks but can be excluded if wanted (did not work for me in the .ipynb files).')

st.write('2. The following code works in theory for any positive N of type int. But in practice, the code takes to long for N>3, so it works fine for one, two and three Pendulums but not for four or more. There is no error in the code but simply a bottleneck: in one step the Lagrange Equations need to be rearanged for the accelerations which - with rising N - causes an explosion in the computational complexity, since every Langrange Equation seems to contain every acceleration (at least observed until N=3) so the function "solve" in the code from sympy is simply not able to rearange the equations: I also tried, as one can see, the "LagrangeMethod" provided by "sympy.physics.mechanics" and try to bypass this bottleneck by using a matrix representation ( and thus trying to invert a matrix) but even this method worked only up N=3 so the Lagrange Methode simply could not handle the problem as I wanted to.')

st.write('3. The following code shows the complete progress from the definition of the parameters until the animation of the numerical results. If one is just interested in the symbolic equations and entities, one should modify the code and only keep the part until the line "time_6 = time.time()". After that the preparation for the numerical simulation starts. One step missing would be calling the function "write_latex_to_txt" to write the calculated expressions to a .txt file. This should greatly reduce runtime and should enable the code to calculate the Lagrange Equations and the other entities for N>3. The follwoing code refers to "N-fold_Pendulum_Lagrange.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Lagrangian, Lagrange Equations, Mass Matrix, Force Vector, Potential Energy and Kinetic Energy')

st.write('(Note: Lagrangian and Lagrange Equations are saved as LaTex expressions to .txt file (see code))')

st.write('(Note: the following Lagrangian, Lagrange Equations, Mass Matrix and Force Vector belong to a system with three pendulums (N=3); for bigger N the code did not finish (see text above.))')

st.write('Lagrangian is:')

st.latex(r'''- g m_{1} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1}\right) - g m_{2} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1} - l_{2} \cos{\left(\phi_{2}{\left(t \right)} \right)} + l_{2}\right) - g m_{3} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1} - l_{2} \cos{\left(\phi_{2}{\left(t \right)} \right)} + l_{2} - l_{3} \cos{\left(\phi_{3}{\left(t \right)} \right)} + l_{3}\right) + \frac{l_{1}^{2} m_{1} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + \frac{l_{1}^{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + \frac{l_{1}^{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + \frac{l_{2}^{2} m_{2} \left(\frac{d}{d t} \phi_{2}{\left(t \right)}\right)^{2}}{2} + \frac{l_{2}^{2} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)}\right)^{2}}{2} + l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + \frac{l_{3}^{2} m_{3} \left(\frac{d}{d t} \phi_{3}{\left(t \right)}\right)^{2}}{2}''')

st.write('Lagrange Equations are:')

st.latex(r'''\left[\begin{matrix}g l_{1} m_{1} \sin{\left(\phi_{1}{\left(t \right)} \right)} + g l_{1} m_{2} \sin{\left(\phi_{1}{\left(t \right)} \right)} + g l_{1} m_{3} \sin{\left(\phi_{1}{\left(t \right)} \right)} + l_{1}^{2} m_{1} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} + l_{1}^{2} m_{2} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} + l_{1}^{2} m_{3} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} - l_{1} l_{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{2}{\left(t \right)} - l_{1} l_{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{2}{\left(t \right)} - l_{1} l_{3} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{3}{\left(t \right)}\\g l_{2} m_{2} \sin{\left(\phi_{2}{\left(t \right)} \right)} + g l_{2} m_{3} \sin{\left(\phi_{2}{\left(t \right)} \right)} - l_{1} l_{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} - l_{1} l_{2} m_{2} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} - l_{1} l_{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} - l_{1} l_{2} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} + l_{2}^{2} m_{2} \frac{d^{2}}{d t^{2}} \phi_{2}{\left(t \right)} + l_{2}^{2} m_{3} \frac{d^{2}}{d t^{2}} \phi_{2}{\left(t \right)} - l_{2} l_{3} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{3}{\left(t \right)}\\g l_{3} m_{3} \sin{\left(\phi_{3}{\left(t \right)} \right)} - l_{1} l_{3} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} - l_{1} l_{3} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{1}{\left(t \right)} - l_{2} l_{3} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} - l_{2} l_{3} m_{3} \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \phi_{2}{\left(t \right)} + l_{3}^{2} m_{3} \frac{d^{2}}{d t^{2}} \phi_{3}{\left(t \right)}\end{matrix}\right]''')

st.write('Mass Matrix is:')

st.latex(r'''\left[\begin{matrix}l_{1}^{2} m_{1} + l_{1}^{2} m_{2} + l_{1}^{2} m_{3} & l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} & l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)}\\l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} & l_{2}^{2} m_{2} + l_{2}^{2} m_{3} & l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)}\\l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} & l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} & l_{3}^{2} m_{3}\end{matrix}\right]''')

st.write('Force Vector is:')

st.latex(r'''\left[\begin{matrix}- g l_{1} m_{1} \sin{\left(\phi_{1}{\left(t \right)} \right)} - g l_{1} m_{2} \sin{\left(\phi_{1}{\left(t \right)} \right)} - g l_{1} m_{3} \sin{\left(\phi_{1}{\left(t \right)} \right)} + l_{1} l_{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} - l_{1} l_{2} m_{2} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} - l_{1} l_{2} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{3} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} - l_{1} l_{3} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)}\\- g l_{2} m_{2} \sin{\left(\phi_{2}{\left(t \right)} \right)} - g l_{2} m_{3} \sin{\left(\phi_{2}{\left(t \right)} \right)} + l_{1} l_{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} + l_{1} l_{2} m_{2} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{2}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} + l_{1} l_{2} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{2} l_{3} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} - l_{2} l_{3} m_{3} \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)}\\- g l_{3} m_{3} \sin{\left(\phi_{3}{\left(t \right)} \right)} + l_{1} l_{3} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} + l_{1} l_{3} m_{3} \sin{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)} - \frac{d}{d t} \phi_{3}{\left(t \right)}\right) \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{2} l_{3} m_{3} \sin{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)}\end{matrix}\right]''')

st.write('Potential Energy U is:')

st.latex(r'''g m_{1} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1}\right) + g m_{2} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1} - l_{2} \cos{\left(\phi_{2}{\left(t \right)} \right)} + l_{2}\right) + g m_{3} \left(- l_{1} \cos{\left(\phi_{1}{\left(t \right)} \right)} + l_{1} - l_{2} \cos{\left(\phi_{2}{\left(t \right)} \right)} + l_{2} - l_{3} \cos{\left(\phi_{3}{\left(t \right)} \right)} + l_{3}\right)''')

st.write('Kinetic Energy T is:')

st.latex(r'''\frac{l_{1}^{2} m_{1} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + \frac{l_{1}^{2} m_{2} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + \frac{l_{1}^{2} m_{3} \left(\frac{d}{d t} \phi_{1}{\left(t \right)}\right)^{2}}{2} + l_{1} l_{2} m_{2} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{2}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\phi_{1}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{1}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + \frac{l_{2}^{2} m_{2} \left(\frac{d}{d t} \phi_{2}{\left(t \right)}\right)^{2}}{2} + \frac{l_{2}^{2} m_{3} \left(\frac{d}{d t} \phi_{2}{\left(t \right)}\right)^{2}}{2} + l_{2} l_{3} m_{3} \cos{\left(\phi_{2}{\left(t \right)} - \phi_{3}{\left(t \right)} \right)} \frac{d}{d t} \phi_{2}{\left(t \right)} \frac{d}{d t} \phi_{3}{\left(t \right)} + \frac{l_{3}^{2} m_{3} \left(\frac{d}{d t} \phi_{3}{\left(t \right)}\right)^{2}}{2}''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code('''
# Code for N-fold pendulum; works only up to N=3
# afterwards solving equations for accelerations takes too long
# no solution was obtained because simulation was aborted after some time; it theoretically works for N>3 but can not be confirmed

# ------------------------------------------------------------------------------------------------------------------------ #

# import all necessary packages
import numpy as np
import sympy as sp

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import matplotlib

import sys
import traceback

from scipy.integrate import solve_ivp

import time

from sympy.physics.mechanics import *

# ------------------------------------------------------------------------------------------------------------------------ #

# in original code functions are defined here; use code from github repository or copy this code and insert them
    
time_1 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# most important parameter: N defines the number of pendulums used
# therefore parameter arrays for m_par, l_par, phi_par and phid_par must be of dimension (1,N)
N = 3

# define variable name here; everything possible; string needed
variable_name = "phi"

# define necessary parameters; m in kg; l in m; phi in degrees; phid in degrees/s; lengths of array must match N
m_par = [0.5, 0.5, 0.5]
l_par = [1.0, 1.0, 1.0]
phi_par = [135, 155, 175]
phid_par = [0, 0, 0]

# ------------------------------------------------------------------------------------------------------------------------ #

# value of g; do not change !!!
g_par = 9.81

# t_start: start simulation time; t_end: end simluation time; steps: number of steps between t_start and t_end
# time vector for simulation and animation; number of samples should not be changed to ensure real-time behaviour of animation and saved gif
# 50 fps chosen as gif could not be displayed with higher values and an increase would not be noticable for humans anyway
t_start = 0
t_end = 10
steps = (t_end-t_start)*50

# ------------------------------------------------------------------------------------------------------------------------ #

# define here method, relative and absolute tolerance of the solve_ivp function used for numerical integration of ODEs
method_ivp = "RK45"
rtol_ivp = 1E-5
atol_ivp = 1E-6

# ------------------------------------------------------------------------------------------------------------------------ #

time_2 = time.time()

# check correct input format of variable_name
check_variable_name(variable_name)
print("")

# check if length parameter arrays matches N
check_length_array(m_par, N)
check_length_array(l_par, N)
check_length_array(phi_par, N)
check_length_array(phid_par, N)
print("")

# check right input format of parameters
check_array_format(m_par, N)
print("")
check_array_format(l_par, N)
print("")
check_array_format(phi_par, N)
print("")
check_array_format(phid_par, N)
print("")

# do necessary precalculations; calculation from degree or degree/s into rad or rad/s
phi_par = array_angular_conversion(phi_par, N)
print("")
phid_par = array_angular_conversion(phid_par, N)
print("")

time_3 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create symbolic parameters
params_m = create_params(m_par, N)
print("Symbolic parameters created for masses are:", params_m)
params_l = create_params(l_par, N)
print("Symbolic parameters created for lengths are:", params_l)
g = [sp.symbols("g")]
print("Symbolic parameter created for standard gravity is:", g)
print("")

# create smybolic time
t = sp.Symbol("t")

# create symbolic variables and time derivatives
variables, derivatives = create_variables_and_derivatives(variable_name, N)
print("Symbolic variables created (in dependence of time) are:", variables)
print("")
print("Symbolic derivatives calculated are:", derivatives)
print("")

time_4 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create additional variables for x and y; x and y are substituted with angles; easier notation of U and especially T in x and y coordinates
x_variables, y_variables = create_x_and_y(N)
print("Symbolic variables created for x are:", x_variables)
print("Symbolic variables created for y are:", y_variables)
print("")

# calculate all xn and yn for every end point of each pendulum; used in next step for T and U notation
xn_variables, yn_variables = create_xn_and_yn(params_l, x_variables, y_variables, N)
print("Symbolic variables created for xn are:", xn_variables)
print("Symbolic variables created for yn are:", yn_variables)
print("")

# kinetic energy T
T = create_T(params_m, g, xn_variables, yn_variables, N)

# simplify T because squaring of sum terms otherwise expands the equation unnecessary
T = T.expand()
T = sp.trigsimp(T)
print("Kinetic energy T was simplified:", T)
print("")

# potential energy
U = create_U(params_m, params_l, g, yn_variables, N)

# determine Lagrange Function L from T and U
L = T - U 

# removes unnecessary factors and digits (in this case factors of 1.0)
L = sp.nsimplify(L)
T = sp.nsimplify(T)
U = sp.nsimplify(U)

# create latex expression of L
L_latex = sp.latex(L)
T_latex = sp.latex(T)
U_latex = sp.latex(U)

print("Lagrangian is:", L)
print("")

time_5 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# LagrangesMethod used to calculate Lagrange Equation, mass matrix and force vector from matrix notation of equations
# equations from manual implementation in next step are used in the rest of the code but for automation for
# other systems these equations should be choosen
LM = LagrangesMethod(L, variables)

# extract equations, mass matrix and force vector
LM_Eq = LM.form_lagranges_equations()
LM_M = LM.mass_matrix
LM_F = LM.forcing

# create latex expression of mass matrix, force vector and Lagrange Equations
LM_M_latex = sp.latex(LM_M)
LM_F_latex = sp.latex(LM_F)
LM_Eq_latex = sp.latex(LM_Eq)

# print mass matrix
print("Mass matrix is:")
print(LM_M)
print("")

# print force vector
print("Force vector is:")
print(LM_F)
print("")

# print lagrange equations
print("Lagrange Equations are:")
print(LM_Eq_latex)
print("")

# replace some wrong symbols in latex expression; otherwise equations and other properties would not show greek letters
L_latex, U_latex, T_latex, LM_M_latex, LM_F_latex, LM_Eq_latex = replace_variable_name(variable_name, L_latex, U_latex, T_latex, LM_M_latex, LM_F_latex, LM_Eq_latex, N)

# ------------------------------------------------------------------------------------------------------------------------ #

# create auxilliary expressions; calculate N Lagrange Equations; display them
Eq = calculate_Lagrange_Equations(L, variable_name, t, N)
print("")

# simplify equations
Eq = simplify_Eq(Eq, N)

time_6 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# list of accelerations
acc = create_acceleration_list(derivatives, N)
print("Accelerations are:", acc)
print("")

time_before_sim = time.time()

# solve equations for acceleration symbols
res = sp.solve(Eq, acc)

time_after_sim = time.time()

print("Solving of equations sucessfully. Time needed:", time_after_sim-time_before_sim, "seconds")
print("")
print("Dimension of res is:", len(res))
print("")

# create phidd expresions
ddt_expressions = create_ddt_expr(res, variable_name, N)
print("")

time_7 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create symbols in tuples without t dependence for variables and derivatives
rplmts = create_symbols_without_t(variable_name, variables, derivatives, N)
print("")

# create symbols in tuples for parameters with their respective numeric value
params_values = create_params_symbols(params_m, params_l, g, m_par, l_par, g_par)
print("")
      
# creation of python functions
func_array = create_python_functions(rplmts, params_values, variable_name, ddt_expressions, N)

time_8 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

tt = np.linspace(t_start, t_end, steps)

# create array start conditions
zz0 = create_zz0(phi_par, phid_par, N)
print("Array zz0 with starting conditions:", zz0)
print("")

# do the numerical integration
res = solve_ivp(rhs, (tt[0], tt[-1]), zz0, t_eval=tt, method=method_ivp, rtol=rtol_ivp, atol=atol_ivp)
print("Numerical integration was sucessful")
print("")

time_9 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# calculate T and U with the rsults from numerical integration
T_sol, U_sol = calculate_T_and_U(variable_name, res, variables, derivatives, T, U, params_values, N)

# Unpacking of individual state components.
var_deg_array = extract_and_convert_res(res, N)

# calculate x and y from solved angle
x_array, y_array = calculate_x_and_y_from_angle(res, l_par, N)

time_10 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create txt file and write tex expressions in it
write_latex_to_txt(L_latex, U_latex, T_latex, LM_Eq_latex, LM_M_latex, LM_F_latex)

time_11 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

print("")
print("Time needed overall is:", time_11 - time_1, "seconds")
print("")
print("Time needed for parameters and function definitions:", time_2 - time_1, "seconds")
print("Time needed for checking of parameters", time_3 - time_2, "seconds")
print("Time needed for creation of symbolic parameters and variables:", time_4 - time_3, "seconds")
print("Time needed for creation of symbolic x and y plus formation of U, T and L:", time_5 - time_4, "seconds")
print("Time needed for calculation of Lagrange Equations, Mass Matrix and Force Vector:", time_6 - time_5, "seconds")
print("Time needed for solving of Equations for accelerations:", time_7 - time_6, "seconds")
print("Time needed for creation of time independent varaibles:", time_8 - time_7, "seconds")
print("Time needed for numerical integration:", time_9 - time_8, "seconds")
print("Time needed for variable conversion:", time_10 - time_9, "seconds")
print("Time needed for saving of variables to txt-file:", time_11 - time_10, "seconds")
print("")

# ------------------------------------------------------------------------------------------------------------------------ #
''')

st.subheader('Output')

st.code('''variable_name is in the right format of str

Array length of m_par ( 3 ) matches N ( 3 )
Array length of l_par ( 3 ) matches N ( 3 )
Array length of phi_par ( 3 ) matches N ( 3 )
Array length of phid_par ( 3 ) matches N ( 3 )

Parameter m_par = 0.5 is in the right format of float
Parameter m_par = 0.5 is in the right format of float
Parameter m_par = 0.5 is in the right format of float

Parameter l_par = 1.0 is in the right format of float
Parameter l_par = 1.0 is in the right format of float
Parameter l_par = 1.0 is in the right format of float

Parameter phi_par = 135 is in the right format of int
Parameter phi_par = 155 is in the right format of int
Parameter phi_par = 175 is in the right format of int

Parameter phid_par = 0 is in the right format of int
Parameter phid_par = 0 is in the right format of int
Parameter phid_par = 0 is in the right format of int

Parameter phi_par = 135 in degrees was transformed to 2.356194490192345 in rad
Parameter phi_par = 155 in degrees was transformed to 2.705260340591211 in rad
Parameter phi_par = 175 in degrees was transformed to 3.0543261909900767 in rad

Parameter phid_par = 0 in degrees was transformed to 0.0 in rad
Parameter phid_par = 0 in degrees was transformed to 0.0 in rad
Parameter phid_par = 0 in degrees was transformed to 0.0 in rad

Symbolic parameters created for masses are: [m1, m2, m3]
Symbolic parameters created for lengths are: [l1, l2, l3]
Symbolic parameter created for standard gravity is: [g]

Symbolic variables created (in dependence of time) are: [phit1(t), phit2(t), phit3(t)]

Symbolic derivatives calculated are: [Derivative(phit1(t), t), Derivative(phit1(t), (t, 2)), Derivative(phit2(t), t), Derivative(phit2(t), (t, 2)), Derivative(phit3(t), t), Derivative(phit3(t), (t, 2))]

Symbolic variables created for x are: [x1(t), x2(t), x3(t)]
Symbolic variables created for y are: [y1(t), y2(t), y3(t)]

Symbolic variables created for xn are: [l1*sin(phit1(t)), l1*sin(phit1(t)) + l2*sin(phit2(t)), l1*sin(phit1(t)) + l2*sin(phit2(t)) + l3*sin(phit3(t))]
Symbolic variables created for yn are: [l1*cos(phit1(t)), l1*cos(phit1(t)) + l2*cos(phit2(t)), l1*cos(phit1(t)) + l2*cos(phit2(t)) + l3*cos(phit3(t))]

Kinetic energy T is: 0.5*m1*(l1**2*sin(phit1(t))**2*Derivative(phit1(t), t)**2 + l1**2*cos(phit1(t))**2*Derivative(phit1(t), t)**2) + 0.5*m2*((-l1*sin(phit1(t))*Derivative(phit1(t), t) - l2*sin(phit2(t))*Derivative(phit2(t), t))**2 + (l1*cos(phit1(t))*Derivative(phit1(t), t) + l2*cos(phit2(t))*Derivative(phit2(t), t))**2) + 0.5*m3*((-l1*sin(phit1(t))*Derivative(phit1(t), t) - l2*sin(phit2(t))*Derivative(phit2(t), t) - l3*sin(phit3(t))*Derivative(phit3(t), t))**2 + (l1*cos(phit1(t))*Derivative(phit1(t), t) + l2*cos(phit2(t))*Derivative(phit2(t), t) + l3*cos(phit3(t))*Derivative(phit3(t), t))**2)

Kinetic energy T was simplified: 0.5*l1**2*m1*Derivative(phit1(t), t)**2 + 0.5*l1**2*m2*Derivative(phit1(t), t)**2 + 0.5*l1**2*m3*Derivative(phit1(t), t)**2 + 1.0*l1*l2*m2*cos(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + 1.0*l1*l2*m3*cos(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + 1.0*l1*l3*m3*cos(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t) + 0.5*l2**2*m2*Derivative(phit2(t), t)**2 + 0.5*l2**2*m3*Derivative(phit2(t), t)**2 + 1.0*l2*l3*m3*cos(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t) + 0.5*l3**2*m3*Derivative(phit3(t), t)**2

Potential energy U is: g*m1*(-l1*cos(phit1(t)) + l1) + g*m2*(-l1*cos(phit1(t)) + l1 - l2*cos(phit2(t)) + l2) + g*m3*(-l1*cos(phit1(t)) + l1 - l2*cos(phit2(t)) + l2 - l3*cos(phit3(t)) + l3)

Lagrangian is: -g*m1*(-l1*cos(phit1(t)) + l1) - g*m2*(-l1*cos(phit1(t)) + l1 - l2*cos(phit2(t)) + l2) - g*m3*(-l1*cos(phit1(t)) + l1 - l2*cos(phit2(t)) + l2 - l3*cos(phit3(t)) + l3) + l1**2*m1*Derivative(phit1(t), t)**2/2 + l1**2*m2*Derivative(phit1(t), t)**2/2 + l1**2*m3*Derivative(phit1(t), t)**2/2 + l1*l2*m2*cos(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m3*cos(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l3*m3*cos(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t) + l2**2*m2*Derivative(phit2(t), t)**2/2 + l2**2*m3*Derivative(phit2(t), t)**2/2 + l2*l3*m3*cos(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t) + l3**2*m3*Derivative(phit3(t), t)**2/2

Mass matrix is:
Matrix([[l1**2*m1 + l1**2*m2 + l1**2*m3, l1*l2*m2*cos(phit1(t) - phit2(t)) + l1*l2*m3*cos(phit1(t) - phit2(t)), l1*l3*m3*cos(phit1(t) - phit3(t))], [l1*l2*m2*cos(phit1(t) - phit2(t)) + l1*l2*m3*cos(phit1(t) - phit2(t)), l2**2*m2 + l2**2*m3, l2*l3*m3*cos(phit2(t) - phit3(t))], [l1*l3*m3*cos(phit1(t) - phit3(t)), l2*l3*m3*cos(phit2(t) - phit3(t)), l3**2*m3]])

Force vector is:
Matrix([[-g*l1*m1*sin(phit1(t)) - g*l1*m2*sin(phit1(t)) - g*l1*m3*sin(phit1(t)) + l1*l2*m2*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit2(t), t) - l1*l2*m2*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m3*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit2(t), t) - l1*l2*m3*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l3*m3*(Derivative(phit1(t), t) - Derivative(phit3(t), t))*sin(phit1(t) - phit3(t))*Derivative(phit3(t), t) - l1*l3*m3*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t)], [-g*l2*m2*sin(phit2(t)) - g*l2*m3*sin(phit2(t)) + l1*l2*m2*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t) + l1*l2*m2*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m3*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t) + l1*l2*m3*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l2*l3*m3*(Derivative(phit2(t), t) - Derivative(phit3(t), t))*sin(phit2(t) - phit3(t))*Derivative(phit3(t), t) - l2*l3*m3*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t)], [-g*l3*m3*sin(phit3(t)) + l1*l3*m3*(Derivative(phit1(t), t) - Derivative(phit3(t), t))*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t) + l1*l3*m3*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t) + l2*l3*m3*(Derivative(phit2(t), t) - Derivative(phit3(t), t))*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t) + l2*l3*m3*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t)]])

Lagrange Equations are:
\left[\begin{matrix}g l_{1} m_{1} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} \right)} + g l_{1} m_{2} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} \right)} + g l_{1} m_{3} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} \right)} + l_{1}^{2} m_{1} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} + l_{1}^{2} m_{2} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} + l_{1}^{2} m_{3} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} - l_{1} l_{2} m_{2} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{2}{\left(t \right)} - l_{1} l_{2} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{2}{\left(t \right)} - l_{1} l_{3} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{3}{\left(t \right)}\\g l_{2} m_{2} \sin{\left(\operatorname{phit}_{2}{\left(t \right)} \right)} + g l_{2} m_{3} \sin{\left(\operatorname{phit}_{2}{\left(t \right)} \right)} - l_{1} l_{2} m_{2} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - l_{1} l_{2} m_{2} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{2} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} - l_{1} l_{2} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - l_{1} l_{2} m_{3} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} + l_{1} l_{2} m_{3} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{2}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} + l_{2}^{2} m_{2} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{2}{\left(t \right)} + l_{2}^{2} m_{3} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{2}{\left(t \right)} - l_{2} l_{3} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \sin{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \cos{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{3}{\left(t \right)}\\g l_{3} m_{3} \sin{\left(\operatorname{phit}_{3}{\left(t \right)} \right)} - l_{1} l_{3} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} - l_{1} l_{3} m_{3} \sin{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{1}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{1} l_{3} m_{3} \cos{\left(\operatorname{phit}_{1}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{1}{\left(t \right)} - l_{2} l_{3} m_{3} \left(\frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} - \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)}\right) \sin{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} - l_{2} l_{3} m_{3} \sin{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d}{d t} \operatorname{phit}_{2}{\left(t \right)} \frac{d}{d t} \operatorname{phit}_{3}{\left(t \right)} + l_{2} l_{3} m_{3} \cos{\left(\operatorname{phit}_{2}{\left(t \right)} - \operatorname{phit}_{3}{\left(t \right)} \right)} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{2}{\left(t \right)} + l_{3}^{2} m_{3} \frac{d^{2}}{d t^{2}} \operatorname{phit}_{3}{\left(t \right)}\end{matrix}\right]

Lagrange Equation 1 was calculated sucessfully and saved:
g*l1*m1*sin(phit1(t)) + g*l1*m2*sin(phit1(t)) + g*l1*m3*sin(phit1(t)) + l1**2*m1*Derivative(phit1(t), (t, 2)) + l1**2*m2*Derivative(phit1(t), (t, 2)) + l1**2*m3*Derivative(phit1(t), (t, 2)) - l1*l2*m2*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit2(t), t) + l1*l2*m2*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m2*cos(phit1(t) - phit2(t))*Derivative(phit2(t), (t, 2)) - l1*l2*m3*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit2(t), t) + l1*l2*m3*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m3*cos(phit1(t) - phit2(t))*Derivative(phit2(t), (t, 2)) - l1*l3*m3*(Derivative(phit1(t), t) - Derivative(phit3(t), t))*sin(phit1(t) - phit3(t))*Derivative(phit3(t), t) + l1*l3*m3*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t) + l1*l3*m3*cos(phit1(t) - phit3(t))*Derivative(phit3(t), (t, 2))

Lagrange Equation 2 was calculated sucessfully and saved:
g*l2*m2*sin(phit2(t)) + g*l2*m3*sin(phit2(t)) - l1*l2*m2*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t) - l1*l2*m2*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m2*cos(phit1(t) - phit2(t))*Derivative(phit1(t), (t, 2)) - l1*l2*m3*(Derivative(phit1(t), t) - Derivative(phit2(t), t))*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t) - l1*l2*m3*sin(phit1(t) - phit2(t))*Derivative(phit1(t), t)*Derivative(phit2(t), t) + l1*l2*m3*cos(phit1(t) - phit2(t))*Derivative(phit1(t), (t, 2)) + l2**2*m2*Derivative(phit2(t), (t, 2)) + l2**2*m3*Derivative(phit2(t), (t, 2)) - l2*l3*m3*(Derivative(phit2(t), t) - Derivative(phit3(t), t))*sin(phit2(t) - phit3(t))*Derivative(phit3(t), t) + l2*l3*m3*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t) + l2*l3*m3*cos(phit2(t) - phit3(t))*Derivative(phit3(t), (t, 2))

Lagrange Equation 3 was calculated sucessfully and saved:
g*l3*m3*sin(phit3(t)) - l1*l3*m3*(Derivative(phit1(t), t) - Derivative(phit3(t), t))*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t) - l1*l3*m3*sin(phit1(t) - phit3(t))*Derivative(phit1(t), t)*Derivative(phit3(t), t) + l1*l3*m3*cos(phit1(t) - phit3(t))*Derivative(phit1(t), (t, 2)) - l2*l3*m3*(Derivative(phit2(t), t) - Derivative(phit3(t), t))*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t) - l2*l3*m3*sin(phit2(t) - phit3(t))*Derivative(phit2(t), t)*Derivative(phit3(t), t) + l2*l3*m3*cos(phit2(t) - phit3(t))*Derivative(phit2(t), (t, 2)) + l3**2*m3*Derivative(phit3(t), (t, 2))


Equations sucessfully simplified

Accelerations are: [Derivative(phit1(t), (t, 2)), Derivative(phit2(t), (t, 2)), Derivative(phit3(t), (t, 2))]

Solving of equations sucessfully. Time needed: 27.029511213302612 seconds

Dimension of res is: 3

Expressions for ddt terms were created


Replacement symbols for variables were created: [(Derivative(phit1(t), (t, 2)), phidd1), (Derivative(phit1(t), t), phid1), (phit1(t), phi1), (Derivative(phit2(t), (t, 2)), phidd2), (Derivative(phit2(t), t), phid2), (phit2(t), phi2), (Derivative(phit3(t), (t, 2)), phidd3), (Derivative(phit3(t), t), phid3), (phit3(t), phi3)]

Replacement symbols for parameters were created: [(m1, 0.5), (l1, 1.0), (m2, 0.5), (l2, 1.0), (m3, 0.5), (l3, 1.0), (g, 9.81)]

Array zz0 with starting conditions: [2.35619449 0.         2.70526034 0.         3.05432619 0.        ]

Numerical integration was sucessful

Unpacking and conversion of simulation results sucessful; order of entries: var1, d/dt var1, ...


Time needed overall is: 35.471991777420044 seconds

Time needed for parameters and function definitions: 0.0015027523040771484 seconds
Time needed for checking of parameters 0.0010480880737304688 seconds
Time needed for creation of symbolic parameters and variables: 0.0031585693359375 seconds
Time needed for creation of symbolic x and y plus formation of U, T and L: 2.499464988708496 seconds
Time needed for calculation of Lagrange Equations, Mass Matrix and Force Vector: 3.537053108215332 seconds
Time needed for solving of Equations for accelerations: 27.029511213302612 seconds
Time needed for creation of time independent varaibles: 1.468013048171997 seconds
Time needed for numerical integration: 0.9307358264923096 seconds
Time needed for variable conversion: 0.0010008811950683594 seconds
Time needed for saving of variables to txt-file: 0.0005033016204833984 seconds''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code('''
for i in range(N):
    # creation of labels
    var_name = variable_name + str(i+1)
    vard_name = variable_name + "d" + str(i+1)
    y_label = var_name +", " + vard_name

    # note for users because automtic labeling with tex symbols did not work
    print("")
    print("Note: " + vard_name + " represents time derivative of " + var_name)
    print("")
  
    # plot of angle and angular velocity over time
    plt.figure(i+1)
    plt.plot(tt, var_deg_array[2*i], 'r', lw=2, label=var_name)
    plt.plot(tt, var_deg_array[2*i+1], 'b', lw=2, label=vard_name)
    plt.title(str(N) + "-fold Pendulum: Angle and Angular Velocity over Time")
    plt.legend()
    plt.xlabel('Time (seconds)')
    plt.ylabel(y_label)
    plt.grid()
    plt.show()

# ------------------------------------------------------------------------------------------------------------------------ #

# plot of kinetic energy T and potential energy U over time
plt.figure(N+1)
plt.plot(tt, T_sol, 'r', lw=2, label='T')
plt.plot(tt, U_sol, 'b', lw=2, label='U')
plt.title('Kinetic and Potential Energy over Time')
plt.legend()
plt.xlabel('Time (seconds)')
plt.ylabel(r'T (Joule), U (Joule)')
plt.grid()
plt.show()

# ------------------------------------------------------------------------------------------------------------------------ #

# calculate mean from sum of T and U -> controll if numerical simulation was sucessful (great shifts would be a hint for bad simulation results)
# mean substracted in following plot to visualize only the deviation from zero -> ideally should be zero but small deviations are okay
mean_T_plus_U = np.mean(T_sol+U_sol)

# plot of deviation of sum of T and U over time
plt.figure(N+2)
plt.plot(tt, T_sol + U_sol - mean_T_plus_U, 'r', lw=2, label='sum of T and U')
plt.title('Sum of kinetic and potential energy over time (minus mean)')
plt.legend()
plt.xlabel('Time (seconds)')
plt.ylabel('sum of T and U')
plt.grid()
plt.show()

# ------------------------------------------------------------------------------------------------------------------------ #

# determine maximum positive and negative deviation; return both values as a positive percent value
max_error = (abs(max(T_sol + U_sol - mean_T_plus_U))/mean_T_plus_U)*100
min_error = (abs(min(T_sol + U_sol - mean_T_plus_U))/mean_T_plus_U)*100

print("")
print("Mean of sum of T and U is:", mean_T_plus_U)
print("Deviation is a maximum of", max_error, "percent in positive direction")
print("Deviation is a maximum of", min_error, "percent in negative direction")
print("")

# ------------------------------------------------------------------------------------------------------------------------ #

# y over x values
plt.figure(N+3)
for i in range(N):
    label_loop = 'x' + str(i+1) + ', y' + str(i+1)
    plt.plot(x_array[i], y_array[i], lw=2, label=label_loop)
title = str(N) + '-fold Pendulum: y over x'
plt.title(title)
plt.legend()
plt.xlabel('x (meter)')
plt.ylabel('y (meter)')
plt.grid()
plt.show()
''')

st.subheader('Output')

st.code('''Note: thetad1 represents time derivative of theta1''')

st.image('page_1_pic_1.png')

st.code('''Note: thetad2 represents time derivative of theta2''')

st.image('page_1_pic_2.png')

st.code('''Note: thetad3 represents time derivative of theta3''')

st.image('page_1_pic_3.png')

st.image('page_1_pic_4.png')

st.image('page_1_pic_5.png')

st.code('''Mean of sum of T and U is: 53.61626279260417
Deviation is a maximum of 0.014345397189405192 percent in positive direction
Deviation is a maximum of 0.007883460483628653 percent in negative direction''')

st.image('page_1_pic_6.png')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code('''# animation of angles over time; angels are drawn in "real-time"

# ------------------------------------------------------------------------------------------------------------------------ #

# import necessary packages
from IPython.display import HTML

# ------------------------------------------------------------------------------------------------------------------------ #

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

# ------------------------------------------------------------------------------------------------------------------------ #

# calculate animation
ani_ang = animate_ang(var_deg_array)

# open animation with this line
HTML(ani_ang.to_jshtml())

# save animation as gif to document results; dpi could be varied; do not change fps value: otherwise no real time behaviour!
from matplotlib.animation import FuncAnimation, PillowWriter
ani_ang.save("N-fold_Pendulum_angle.gif", dpi=300, writer=PillowWriter(fps=50))

''')

st.subheader('Output')

col1, col2 = st.columns(2)

with col1:
    st.image('page_1_pic_7.png')
    
with col2:
    import base64
    file = open(r"N-fold_Pendulum_angle.gif", "rb")
    contents = file.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file.close()

    st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
    
# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code('''
# animation of angular velocity over time; angular velocitys are drawn in "real-time"

# ------------------------------------------------------------------------------------------------------------------------ #

# import necessary packages
from IPython.display import HTML

# ------------------------------------------------------------------------------------------------------------------------ #

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

# ------------------------------------------------------------------------------------------------------------------------ #

# calculate animation
ani_ang_vel = animate_ang_vel(var_deg_array)

# open animation with this line
HTML(ani_ang_vel.to_jshtml())

# save animation as gif to document results; dpi could be varied; do not change fps value: otherwise no real time behaviour!
from matplotlib.animation import FuncAnimation, PillowWriter
ani_ang_vel.save("N-fold_Pendulum_angular_velocity.gif", dpi=300, writer=PillowWriter(fps=50))
''')

st.subheader('Output')

col1, col2 = st.columns(2)

with col1:
    st.image('page_1_pic_8.png')
    
with col2:
    import base64
    file = open(r"N-fold_Pendulum_angular_velocity.gif", "rb")
    contents = file.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file.close()

    st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
    
# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code('''
# animation of pendulum motion over time; motion is drawn in "real-time"

# ------------------------------------------------------------------------------------------------------------------------ #

# import necessary packages 
from IPython.display import HTML

# ------------------------------------------------------------------------------------------------------------------------ #

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

# extract angels from res.y and write in array var_array
var_array = []
for i in range(N):
    var_temp = res.y[2*i]
    var_array.append(var_temp)

# calculates sum of all length; used for scaling of plot
sum_length = sum(l_par)

# calculates all x and y values from angles and lengths
x_array, y_array = pend_pos(var_array)

# calculate animation
anim = animate_pendulum(x_array, y_array)

# open animation with this line
HTML(anim.to_jshtml())

# save animation as gif to document results; dpi could be varied; do not change fps value: otherwise no real time behaviour!
from matplotlib.animation import FuncAnimation, PillowWriter
anim.save("N-fold_Pendulum_motion.gif", dpi=300, writer=PillowWriter(fps=50))

''')

st.subheader('Output')

col1, col2 = st.columns(2)

with col1:
    st.image('page_1_pic_9.png')
    
with col2:
    import base64
    file = open(r"N-fold_Pendulum_motion.gif", "rb")
    contents = file.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file.close()

    st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
    

