
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Introduction Kane Methode")

st.write('The following page should present the basic idea of the Kane Methode as developed or presented by Thomas R. Kane.')

st.write('The implementation for this project can be found under sympy.physics.mechanics here: https://docs.sympy.org/latest/modules/physics/mechanics/kane.html')

st.write('The first source for this methode is a book published by Kane and others and can be found in Google Scholar under: Kane, Thomas R., and David A. Levinson. Dynamics, Theory and Applications. New York: McGraw-Hill, 1985. Print.')

st.write('Since the main focus of this project was not the description of the mathematical framework for the Kane Methode I will simply present a brief description in the next section as fas as I understood the concept. Interested readers are referred to the sources mentioned above.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Description Kane Methode')

st.write('Kanes method, named after Thomas R. Kane, is a powerful technique in analytical mechanics for deriving the equations of motion for mechanical systems. This method is especially advantageous for systems with complex constraints and nonholonomic conditions. Kanes method simplifies the process of obtaining equations of motion, making it particularly useful for systems with multiple interconnected parts.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Key Features Kane Methode')

st.write('Kanes method introduces several concepts that distinguish it from traditional approaches like the Lagrangian method. One of the primary differences is the use of generalized speeds instead of generalized coordinates. Generalized speeds are directly related to the systems motion and are not necessarily the time derivatives of generalized coordinates. This approach often makes the resulting equations more intuitive and manageable.')

st.write('To begin with, Kanes method establishes kinematic differential equations that relate the generalized speeds to the generalized coordinates and their time derivatives. A critical component in this process is the concept of partial velocities, which are the velocities of the systems particles with respect to the generalized speeds. Partial velocities simplify the expressions for kinetic energy and generalized forces, thereby streamlining the derivation of the equations of motion.')

st.write('The method results in a set of dynamic equations that are typically fewer and simpler than those obtained through Lagrangian or Newtonian methods. This efficiency is particularly beneficial for complex systems. Furthermore, Kanes method effectively handles nonholonomic systems—systems with constraints that depend on velocities—by directly incorporating these constraints into the formulation of the equations of motion.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Steps in Kane Methode')

st.write('The process of using Kanes method involves several key steps:')

latex_q = 'q_i'
latex_u = 'u_i'

st.write(f'1. Identify Generalized Coordinates and Speeds: Define the generalized coordinates ${latex_q}$ to describe the systems configuration and the generalized speeds ${latex_u}$ to describe its motion.')

st.write('2. Express Velocities: Determine the velocities of all particles in the system in terms of the generalized coordinates and speeds.')

st.write('3. Calculate Partial Velocities: Find the partial velocities of each particle with respect to each generalized speed.')

st.write('4. Formulate Kinematic Equations: Establish the kinematic differential equations relating the time derivatives of the generalized coordinates to the generalized speeds.')

st.write('5. Compute Forces and Moments: Determine the generalized active forces and moments acting on the system.')

st.write('6. Apply Kane’s Equations: Use Kane’s equations to relate the generalized active forces and moments to the partial velocities, resulting in the equations of motion.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Advantages of Kane Methode')

st.write('Kanes method is known for its efficiency, often resulting in fewer and simpler equations than traditional methods. This makes it computationally efficient and well-suited for systems with complex kinematics and constraints. It is particularly effective in fields such as robotics, vehicle dynamics, aerospace engineering, and biomechanics.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('Example: Double Pendulum')

latex_l1 = 'l_1'
latex_l2 = 'l_2'
latex_m1 = 'm_1'
latex_m2 = 'm_2'
latex_theta1 = r'''\theta_1'''
latex_theta2 = r'''\theta_2'''
latex_u1 = 'u_1'
latex_u2 = 'u_2'
latex_thetad1 = r'''\dot{\theta_1}'''
latex_thetad2 = r'''\dot{\theta_2}'''
latex_v1 = 'v_1'
latex_v2 = 'v_2'
latex_vec1 = r'''l_1  \dot{\theta_1}'''
latex_vec2 = r'''l_1  \dot{\theta_1} + l_2  \dot{\theta_2}'''

st.write(f'Consider a double pendulum, which consists of two pendulums attached end to end. The first pendulum has length ${latex_l1}$ and mass ${latex_m1}$ and the second pendulum has length ${latex_l2}$ and mass ${latex_m2}$')

st.write(f'1. Generalized Coordinates and Speeds: Let ${latex_theta1}$ and ${latex_theta2}$ be the angles of the first and second pendulums, respectively, with the vertical. Define the generalized speeds ${latex_u1}$ = ${latex_thetad1}$ and ${latex_u2}$ = ${latex_thetad2}$.')

st.write(f'2. Velocities: The velocity of the first mass ${latex_m1}$ is given by: ${latex_v1}$ = ${latex_vec1}$. The velocity of the second mass ${latex_m2}$ is given by: ${latex_v2}$ = ${latex_vec2}$. Note that ${latex_v1}$ and ${latex_v2}$ are vectorized values.')

st.write(f'3. Partial Velocities: The partial velocities for ${latex_m1}$ and ${latex_m2}$ with respect to ${latex_thetad1}$ and ${latex_thetad2}$ are calculated.')

st.write(f'4. Kinematic Equations: The kinematic differential equations relate ${latex_thetad1}$ and ${latex_thetad2}$ to the generalized speeds ${latex_u1}$ and ${latex_u2}$')

st.write(f'5. Forces and Moments: Consider the gravitational forces acting on ${latex_m1}$ and ${latex_m2}$. These are incorporated into the generalized active forces.')

st.write(f'6. Kanes Equations: Apply Kanes equations to derive the equations of motion. The resulting equations describe the dynamics of the double pendulum in terms of ${latex_theta1}$, ${latex_theta2}$, ${latex_u1}$ and ${latex_u2}$')

st.write('By following these steps, Kanes method provides a systematic approach to deriving the equations of motion for the double pendulum, illustrating its power and efficiency in handling complex mechanical systems.')


