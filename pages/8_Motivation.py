
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Motivation")

st.write('The following text is the motivation letter for my project idea translated into english. It describes the goal, the intermediate targets and the way it should be implemented.')

st.subheader('Motivation Letter')

st.write('As an engineer, the modelling and simulation of physical systems is a recurring part of your daily work. One type of system is a mechanical system, which can be abstracted using components such as rigid bodies or masses, springs and dampers. This allows a simpler mathematical consideration of the behaviour, which can be expressed by equations of motion in dynamic systems. One approach is the Lagrange formalism, which uses the Lagrange equations of the second kind to describe systems by means of their kinetic and potential energy and a series of generalised coordinates.')

st.write('This approach will be applied to pendulums in the following project. In the simplest form, a massless rod is considered, at one end of which a point mass is assumed and the other end is fixed with a suitable bearing. If you want to pursue this approach further, you can consider several coupled pendulums. In the case of the project, a special case of this will be considered: each additional pendulum is attached to the free end of the previous one, whereby the basis for each individual pendulum is assumed as described above (point mass at one end, suitable bearing at the other end). In addition, the Kane method was selected and implemented as an alternative in the project. Specifically, a solution is to be implemented for any number of pendulums N, their parameters (mass and length) and the subsequent graphical evaluation in a suitable form, for example by means of an animation or other diagrams. In addition, the equations of motion and the formulae for kinetic and potential energy are to be calculated as symbolic equations and saved in LaTex format.')

st.write('The basic procedure is carried out in several steps: The starting point is the description of a simple pendulum, which is used as an example to test the implementation of the Lagrange equations and the graphical representation. This is followed by the step to two pendulums, as the equations of motion can still be set up easily by hand and any errors in the implementation can be detected more easily. Finally, as mentioned above, the step to N pendulums takes place, whereby the challenge here is to implement all steps using suitable functions for arbitrary dimensions N. In addition, the chaotic behaviour of the system should be represented graphically. Once the concept works for N pendulums, there is a basis for future work to adapt this program for other systems, as a modular structure of the individual steps is important.')

st.image('Lagrange_1.png')

st.write('Example for the Lagrange Methode for a double pendulum; you can see here the Lagrangian for an already linearized system (https://www.youtube.com/watch?v=Z7gxaC85JxU)')