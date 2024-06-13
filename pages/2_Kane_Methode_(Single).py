
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Kane Methode (Single)")

st.write('The following page shows the code for a N-fold Pendulum calculated with the Kane Methode. In the next section I will show the code and afterwards the results. The code itself is commented but the functions themself left out for better visibility. Documentation for the functions can be found under "Documentation Kane Methode". I will show in an alternating sequence code snippets (as implemented in jupyter notebook) and the corresponding output. The code itself can be found in the github repository linked in the "Homepage".')

st.write('Before the code and the results themself some remarks:')

st.write('1. I excluded the code for the functions in the snippets below for better visibility. The code for the functions is still inserted in the jupyter notebooks but can be excluded if wanted (did not work for me in the .ipynb files).')

st.write('2. The following code works for any positive N of type int. But in practice, the code takes longer for higher numbers of pendulums. In this case (like in every other presented program besides the N-fold Lagrange Methode) the part of the codes that takes the most time is the presentation of the animation. Even the calculation of the animation itself is pretty fast, but displaying it takes the most time for smaller N. Maybe this could speed up if you use this code in an local IDE or do some improvements of the animation part.')

st.write('3. The following code shows the complete progress from the definition of the parameters until the animation of the numerical results. It refers to the file "Single_n-fold_Pendulum_Kane_Method".')

st.write('4. The idea for these functions stems from a website here: https://jakevdp.github.io/blog/2017/03/08/triple-pendulum-chaos/. I used them as a base and modified them according to my needs.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code(r'''
# Code for a single n-fold Pendulum; Kanes Method is used here

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

# n is number of segments of the pendulum
n=3

# either random is 0 or 1
# random=1: choose a random start angle, mass and length for every pendulum; masses and length between 0 and 1; start angle between 0 and 360 degrees
# random=0: one random start angle is chosen that is applied to all pendulum segments; for this case every segment has length of 1m and mass of 1kg
random = 0

# time vector for simulation and animation; number of samples should not be changed to ensure real-time behaviour of animation and saved gif
# 50 fps chosen as gif could not be displayed with higher values and an increase would not be noticable for humans anyway
t_start = 0
t_end = 10
times = np.linspace(t_start, t_end, (t_end-t_start)*50)

# scaling factors for length and mass; change if necessary; if not used set to 1
# usefull because random function often produces values near to zero which is not good for nice simulation results
length_scaling = 2.5
mass_scaling = 5

# set data limit of animation; in MB; might need to be enlarged for other applications
matplotlib.rcParams['animation.embed_limit'] = 100.0

time_1 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# predefine empty lists and parameters
angles = []
velocities = []
lengths = []
masses = []
sum_length = 0

if random==1:
    # as decribed above: assign random angle, mass and length for every pendulum segment; if necessary use scaling factors
    # append values to lists; add sum_length in each cycle up; used for scaling of plot and animation
    for i in range(n):
        rand_angles = rd.random()
        rand_lengths = rd.random() * length_scaling
        rand_masses = rd.random() * mass_scaling
        angles.append(rand_angles*360)
        velocities.append(0)
        lengths.append(rand_lengths)
        masses.append(rand_masses)
        sum_length = sum_length + rand_lengths
else:
    # one random angle for all pendulums
    rand_angles = rd.random()

    #create lists and sum_length as mentioned above
    for i in range(n):
        angles.append(rand_angles*360)
        velocities.append(0)
        lengths.append(1)
        masses.append(1)
        sum_length = sum_length + 1

time_2 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# integrate pendulum and create p
p = integrate_pendulum(n, times, angles, velocities, lengths, masses)

time_3 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# get x and y coordinates from p and lengths
x, y = get_xy_coords(p, lengths)

time_4 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# reform x and y values so that e.g. all x values for one pendulum segment are in one list, which gots appended to x_list; same for y
x_list = []
y_list = []
for o in range(1,n+1):
    x_temp = []
    y_temp = []
    for i in range(len(times)):
        x_temp.append(x[i][o])
        y_temp.append(y[i][o])
    x_list.append(x_temp)
    y_list.append(y_temp)

# plot y over x values; plot each x and y entry seperatly for adding of label
for i in range(n):
    label_loop = "x" + str(i+1) + ", y" + str(i+1)
    plt.plot(x_list[i], y_list[i], lw=2, label=label_loop)

title = str(n) + '-fold Pendulum'
plt.title(title)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.ylim(-sum_length, sum_length)
plt.xlim(-sum_length, sum_length)
plt.grid()
plt.legend()
plt.show()

time_5 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# create animation
anim = animate_pendulum(x, y)

time_6 = time.time()

# ------------------------------------------------------------------------------------------------------------------------ #

# evaluation of times; not necessary for function
times_initialization = time_2 - time_1
times_integration = time_3 - time_2
times_x_y_conversion = time_4 - time_3
times_plotting = time_5 - time_4
times_animation = time_6 - time_5

print("")
print("Time for initialization:", times_initialization, "seconds")
print("Time for integration:", times_integration, "seconds")
print("Time for x y conversion:", times_x_y_conversion, "seconds")
print("Time for plotting:", times_plotting, "seconds")
print("Time for animation calculation:", times_animation, "seconds")
print("")
print("Time overall:", time_6 - time_1)
print("")

print("Angles used:", angles)
print("Velocities used:", velocities)
print("Lengths used:", lengths)
print("Masses used:", masses)
print("")
''')

st.subheader('Output')

st.image('page_2_pic_1.png')

st.code(r'''Time for initialization: 0.0 seconds
Time for integration: 0.6550853252410889 seconds
Time for x y conversion: 0.0 seconds
Time for plotting: 0.15394258499145508 seconds
Time for animation calculation: 0.0645136833190918 seconds

Time overall: 0.8735415935516357

Angles used: [182.74762122623778, 182.74762122623778, 182.74762122623778]
Velocities used: [0, 0, 0]
Lengths used: [1, 1, 1]
Masses used: [1, 1, 1]''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Code')

st.code(r'''
# execute/open animation with this line
HTML(anim.to_jshtml())

# create name for give in dependence of n
name_gif = "Single_" + str(n) + "-fold_Pendulum.gif"

# save animation as gif to document results; dpi could be varied; do not change fps value: otherwise no real time behaviour!
from matplotlib.animation import FuncAnimation, PillowWriter
anim.save(name_gif, dpi=300, writer=PillowWriter(fps=50))
''')

st.subheader('Output')

col1, col2 = st.columns(2)

with col1:
    st.image('page_2_pic_2.png')
    
with col2:
    import base64
    file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_3-fold_Pendulum.gif", "rb")
    contents = file.read()
    data_url = base64.b64encode(contents).decode("utf-8")
    file.close()

    st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Further Examples')

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_2-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_4-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_6-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_8-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)

st.write('')

import base64
file = open(r"C:\Users\jonas\OneDrive\Dokumente\Uni\Master\4. Semester (SS2024)\Python_Kurs\Projekt\Documentation_website\Single_10-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)


