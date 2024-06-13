
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Homepage")

st.write('This website serves as a documentation for the project "Pysaar 2024". All code produced during this project can be found in github: https://github.com/jonasmartin82/PySaar_2024.git')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Project Description')

st.write('The detailed description of the project can be found under "Motivation". In essence the scope of this project was to calculate the equations of motion for a N-fold pendulum. N-fold in this case refers to N pendulums connected in a row. Furthermore an extraction of the equations of motion in LaTex was desired. For solving this task two approaches were chosen: the Lagrange Methode and the Kane Methode. These two methodes (or their implementation) are shown in the seperate pages "Lagrange Methode", "Kane Methode (Single)" and "Kane Methode (Multiple)". Just choose from the side bar on the left. Also the documentation for the used functions could be found there under "Documentation Lagrange Methode" and "Documentation Kane Methode". As a Bonus I added the code for a Single Pendulum and a Double Pendulum. Both programms emerged during the work at the project and can be used to easily understand the Lagrange Methode and adjust them for other mechanical systems. Also added is the page "Motivation" and the page "Introduction Kane Methode" in which I briefly introduce the Kane Methode.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Examples of Results')

st.write('An example for this implementation could be seen below. It shows an 3-fold pendulum and the resulting motion according to the initial conditions. This animation was derived from the Kane Methode. One can notice the layout of the coordinate system: the x axis is horizontal and the y axis vertical while positive values point to the right respectively to the top.')

import base64
file = open(r"Single_3-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
 
st.write('Furthermore the chaotic characteristic of the system was observed and presented in form of an animation. One could see that even small changes in the initial conditions could lead to a totally different behaviour (shown here with the plot below). The starting position of all the Pendulums was in the upper right corner. For further details look into "Kane Methode (Multiple)".')

st.image('Homepage_pic_2.png')

st.write('The animation below shows the moment when the motion of the Pendulums start to differ from one another. From there on it took only moments until the system breaks into total chaos.')

import base64
file = open(r"Multiple_3-fold_Pendulum.gif", "rb")
contents = file.read()
data_url = base64.b64encode(contents).decode("utf-8")
file.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gif" width="600" height="470">', unsafe_allow_html=True,)
    
# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Further Remarks')

st.write('Besides the final implementations in "N-fold_Pendulum_Lagrange", "Single_n-fold_Pendulum_Kane_Method" and "Multiple_n-fold_Pendulum_Kane_Method" some simple examples for a Single Pendulum and Double Pendulum (both hard coded) solved with the Euler Lagrange Equations are added to the github repository.')


