
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Kane Methode (Multiple)")

st.write('The following page shows the code for multiple N-fold Pendulums calculated with the Kane Methode. In the next section I will show the code and afterwards the results. The code itself is commented but the functions themself left out for better visibility. Documentation for the functions can be found under "Documentation Kane Methode". I will show in an alternating sequence code snippets (as implemented in jupyter notebook) and the corresponding output. The code itself can be found in the github repository linked in the "Homepage".')

st.write('Before the code and the results themself some remarks:')

st.write('1. I excluded the code for the functions in the snippets below for better visibility. The code for the functions is still inserted in the jupyter notebooks but can be excluded if wanted (did not work for me in the .ipynb files).')

st.write('2. The following code works for any positive N and n_pendulums of type int. But in practice, the code takes longer for higher numbers of pendulums. In this case (like in every other presented program besides the N-fold Lagrange Methode) the part of the codes that takes the most time is the presentation of the animation. Even the calculation of the animation itself is pretty fast, but displaying it takes the most time for smaller N. Maybe this could speed up if you use this code in an local IDE or do some improvements of the animation part.')

st.write('3. The following code shows the complete progress from the definition of the parameters until the animation of the numerical results. It refers to the file "Multiple_n-fold_Pendulum_Kane_Method".')

st.write('4. The idea for these functions stems from a website here: https://jakevdp.github.io/blog/2017/03/08/triple-pendulum-chaos/. I used them as a base and modified them according to my needs.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code(r'''
# Code for multiple n-fold Pendulums; Kanes Method is used here

# ------------------------------------------------------------------------------------------------------------------------ #

# import all necessary packages
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import time
import matplotlib
import random as rd

from sympy import symbols
from sympy.physics import mechanics

from sympy import Dummy, lambdify
from scipy.integrate import odeint

from matplotlib import animation

from IPython.display import HTML

# ------------------------------------------------------------------------------------------------------------------------ #

# in original code functions are defined here; use code from github repository or copy this code and insert them

# ------------------------------------------------------------------------------------------------------------------------ #

# number n of segments per pendulum
n=3

# number of pendulums
n_pendulums = 5

# error between each starting position; is added e-times to random start angle with e being the index in for loop which iterates over n_pendulums
error = 1E-6

# time vector for simulation and animation; number of samples should not be changed to ensure real-time behaviour of animation and saved gif
# 50 fps chosen as gif could not be displayed with higher values and an increase would not be noticable for humans anyway
t_start = 0
t_end = 10
times = np.linspace(t_start, t_end, (t_end-t_start)*50)

# either predefine rand_angles with desired values in rad or use rd.random() function; comment lines out respectively
#rand_angles = rd.random()
rand_angles = 0.35

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

time_1 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# predefine empty lists
x_array = []
y_array = []

for e in range(n_pendulums):

    # predefine empty lists and parameters
    angles = []
    velocities = []
    lengths = []
    masses = []
    sum_length = 0

    # use the same random start angle for each pendulum and each pendulum segment; same for mass and length; fixed at 1m and 1kg
    # use of fixed values here for simplification of animation; with all variables random it could look to chaotic
    # do display the chaotic nature of the system the little error defined above is enough
    for i in range(n):
        angles.append((rand_angles+e*error)*360)
        velocities.append(0)
        lengths.append(1)
        masses.append(1)
        sum_length = sum_length + 1

    # integrate pendulum and create p
    p = integrate_pendulum_multiple(n, times, angles, velocities, lengths, masses)

    # get x and y coordinates from p and lengths
    x, y = get_xy_coords(p, lengths)

    # append calculated x and y for each pendulum (n_pendulums times) to lists
    x_array.append(x)
    y_array.append(y)

time_2 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create title for plot
title = str(n_pendulums) + ' Pendulums (' + str(n) + '-fold)'

# plot each pair of x and y coordinates
for i in range(n_pendulums):
    plt.plot(x_array[i], y_array[i], lw=2)
   
plt.title(title)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.grid()
plt.show()

time_3 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

# create animation
anim = animate_pendulum_multiple(x_array, y_array)

time_4 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# evaluation of times; not necessary for function
times_simulation = time_2 - time_1
times_plotting = time_3 - time_2
times_animation_calculation = time_4 - time_3

print("")
print("Time for simulation:", times_simulation, "seconds")
print("Time for plotting:", times_plotting, "seconds")
print("Time for animation calculation:", times_animation_calculation, "seconds")
print("")
print("Time overall:", time_4 - time_1, "seconds")
print("")
print("Angle used was:", rand_angles)
print("")''')

st.subheader('Output')

st.image('page_3_pic_1.png')

st.code(r'''Time for simulation: 3.431399345397949 seconds
Time for plotting: 0.2041311264038086 seconds
Time for animation calculation: 0.07972884178161621 seconds

Time overall: 3.715259313583374 seconds

Angle used was: 0.35''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code(r'''
# execute/open animation with this line
HTML(anim.to_jshtml())

# create name for give in dependence of n
name_gif = "Multiple_" + str(n) + "-fold_Pendulum.gif"

# save animation as gif to document results; dpi could be varied; do not change fps value: otherwise no real time behaviour!
from matplotlib.animation import FuncAnimation, PillowWriter
anim.save(name_gif, dpi=300, writer=PillowWriter(fps=50))
''')

st.subheader('Output')

col1, col2 = st.columns(2)

with col1:
    st.image('page_3_pic_2.png')
    
with col2:
    import base64
    file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Multiple_3-fold_Pendulum.gif", "rb")
    contents = file.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file.close()

    st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
    
# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Further Examples')

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Multiple_2-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Multiple_4-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Multiple_6-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)




