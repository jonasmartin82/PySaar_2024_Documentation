
import streamlit as st 

st.set_page_config(page_title = "Python Project Documentation (PySaar 2024)", layout = "wide")

st.logo('logo.png')

st.image('logo.png')

st.title("Documentation Kane Methode")

st.write('In this page the functions used in "Single_n-fold_Pendulum_Kane_Methode" and "Multiple_n-fold_Pendulum_Kane_Mehode" are described. General information about the Kane Methode respectivaly its implementation could be found in the page "Kane Methode (Single)" and "Kane Methode (Multiple)".')

st.write('The order of the functions is like in the jupyter notebook and in the same order as the functions are used. There will be a heading for every function followed by the code with an included docstring as the functions were documented in such a way.')

st.write('Note: the idea for these functions stems from a website here: https://jakevdp.github.io/blog/2017/03/08/triple-pendulum-chaos/. I used them as a base and modified them according to my needs.')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('integrate_pendulum')

st.code('''def integrate_pendulum(n, times, initial_positions, initial_velocities, lengths, masses):
    """
    Uses the Kane Methode to solve the equations of motion for a pendulum with n segments
    Define number of segments, simulation time, initial angles and angular velocities aswell as masses and lengths of each segment
    Specific steps are highlightes below with comments in code itself

    Args:
        n (int): number of segments of pendulum
        times (numpy array): numpy array of simulation time; defined over start and end time and number of steps
        initial_positions (list): list of initial angles of every pendulum segment
        initial_velocities (list): list of initial angular velocities of every pendulum segment
        lengths (list): list of lengths of every pendulum segment
        masses (list): list of masses of every pendulum segment

    Returns:
        p (): points p
    """

    #-------------------------------------------------
    # Step 1: construct the pendulum model
    
    # Generalized coordinates and velocities
    # (in this case, angular positions & velocities of each mass) 
    q = mechanics.dynamicsymbols('q:{0}'.format(n))
    u = mechanics.dynamicsymbols('u:{0}'.format(n))

    # mass and length
    m = symbols('m:{0}'.format(n))
    l = symbols('l:{0}'.format(n))

    # gravity and time symbols
    g, t = symbols('g,t')
    
    #--------------------------------------------------
    # Step 2: build the model using Kane's Method

    # Create pivot point reference frame
    A = mechanics.ReferenceFrame('A')
    P = mechanics.Point('P')
    P.set_vel(A, 0)

    # lists to hold particles, forces, and kinetic ODEs
    # for each pendulum in the chain
    particles = []
    forces = []
    kinetic_odes = []

    for i in range(n):
        # Create a reference frame following the i^th mass
        Ai = A.orientnew('A' + str(i), 'Axis', [q[i], A.z])
        Ai.set_ang_vel(A, u[i] * A.z)

        # Create a point in this reference frame
        Pi = P.locatenew('P' + str(i), l[i] * Ai.x)
        Pi.v2pt_theory(P, A, Ai)

        # Create a new particle of mass m[i] at this point
        Pai = mechanics.Particle('Pa' + str(i), Pi, m[i])
        particles.append(Pai)

        # Set forces & compute kinematic ODE
        forces.append((Pi, m[i] * g * A.x))
        kinetic_odes.append(q[i].diff(t) - u[i])

        P = Pi

    # Generate equations of motion
    KM = mechanics.KanesMethod(A, q_ind=q, u_ind=u, kd_eqs=kinetic_odes)
    fr, fr_star = KM.kanes_equations(particles, forces)
    
    #-----------------------------------------------------
    # Step 3: numerically evaluate equations and integrate

    # initial positions and velocities – assumed to be given in degrees
    y0 = np.deg2rad(np.concatenate([np.broadcast_to(initial_positions, n), np.broadcast_to(initial_velocities, n)]))
        
    # lengths and masses
    lengths = np.broadcast_to(lengths, n)
    masses = np.broadcast_to(masses, n)

    # Fixed parameters: gravitational constant, lengths, and masses
    parameters = [g] + list(l) + list(m)
    parameter_vals = [9.81] + list(lengths) + list(masses)

    # define symbols for unknown parameters
    unknowns = [Dummy() for i in q + u]
    unknown_dict = dict(zip(q + u, unknowns))
    kds = KM.kindiffdict()

    # substitute unknown symbols for qdot terms
    mm_sym = KM.mass_matrix_full.subs(kds).subs(unknown_dict)
    fo_sym = KM.forcing_full.subs(kds).subs(unknown_dict)

    # create functions for numerical calculation 
    mm_func = lambdify(unknowns + parameters, mm_sym)
    fo_func = lambdify(unknowns + parameters, fo_sym)

    # function which computes the derivatives of parameters
    def gradient(y, t, args):
        vals = np.concatenate((y, args))
        sol = np.linalg.solve(mm_func(*vals), fo_func(*vals))
        return np.array(sol).T[0]

    # ODE integration
    return odeint(gradient, y0, times, args=(parameter_vals,))''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('get_xy_coords')

st.code('''def get_xy_coords(p, lengths):
    """
    Get x and y coordinates from lengths of pendulums and p

    Args:
        p (): points p
        lengths (list): list of lengths of each pendulum

    Returns:
        x (): numpy cumsum with all x values
        y (): numpy cumsum with all y values
    """
    # calculate x and y coordinates from lengths and points; safe values in x and y array for animation
    p = np.atleast_2d(p)
    n = p.shape[1] // 2
    lengths_neg = [ -x for x in lengths]
    zeros = np.zeros(p.shape[0])[:, None]
    x = np.hstack([zeros, lengths * np.sin(p[:, :n])])
    y = np.hstack([zeros, lengths_neg * np.cos(p[:, :n])])
    return np.cumsum(x, 1), np.cumsum(y, 1)''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('animate_pendulum')

st.code('''def animate_pendulum(x, y):
    """
    Animates pendulum motion over x and y values; the figure and axes object are defined; creation of line object; definition of init and animate function
    Animate function: for each frame i every line is updated with values for this specific index i

    Args:
        x (List): list of all x values for every pendulum
        y (list): list of all y values for every pendulum

    Returns:
        anim (): animation object
    """
    # define plot and ax; set some parameters
    title = 'Pendel motion over time (Single Pendulum, ' + str(n) + "-fold)"
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')
    ax.set_xlim(-1.05*sum_length, 1.05*sum_length)
    ax.set_ylim(-1.05*sum_length, 1.05*sum_length)
    ax.set_title(title)
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    ax.grid(1)

    # initiate the line object here
    line, = ax.plot([], [], 'o-', lw=2)

    # initialize the line object
    def init():
        line.set_data([], [])
        return line,

    # animate the position of each line for every frame i
    def animate(i):
        line.set_data(x[i], y[i])
        return line,

    # create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(times), interval=1000 * times.max() / len(times), blit=True, init_func=init)

    # closes the static plot which is not used here
    plt.close(fig)
    
    return anim''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('integrate_pendulum_multiple')

st.code('''def integrate_pendulum_multiple(n, times, initial_positions, initial_velocities, lengths, masses):
    """
    Uses the Kane Methode to solve the equations of motion for a pendulum with n segments
    Define number of segments, simulation time, initial angles and angular velocities aswell as masses and lengths of each segment
    Specific steps are highlightes below with comments in code itself

    Args:
        n (int): number of segments of pendulum
        times (numpy array): numpy array of simulation time; defined over start and end time and number of steps
        initial_positions (list): list of initial angles of every pendulum segment
        initial_velocities (list): list of initial angular velocities of every pendulum segment
        lengths (list): list of lengths of every pendulum segment
        masses (list): list of masses of every pendulum segment

    Returns:
        p (): points p
    """
    #-------------------------------------------------
    # Step 1: construct the pendulum model
    
    # Generalized coordinates and velocities
    # (in this case, angular positions & velocities of each mass) 
    q = mechanics.dynamicsymbols('q:{0}'.format(n))
    u = mechanics.dynamicsymbols('u:{0}'.format(n))

    # mass and length
    m = symbols('m:{0}'.format(n))
    l = symbols('l:{0}'.format(n))

    # gravity and time symbols
    g, t = symbols('g,t')
    
    #--------------------------------------------------
    # Step 2: build the model using Kane's Method

    # Create pivot point reference frame
    A = mechanics.ReferenceFrame('A')
    P = mechanics.Point('P')
    P.set_vel(A, 0)

    # lists to hold particles, forces, and kinetic ODEs
    # for each pendulum in the chain
    particles = []
    forces = []
    kinetic_odes = []

    for i in range(n):
        # Create a reference frame following the i^th mass
        Ai = A.orientnew('A' + str(i), 'Axis', [q[i], A.z])
        Ai.set_ang_vel(A, u[i] * A.z)

        # Create a point in this reference frame
        Pi = P.locatenew('P' + str(i), l[i] * Ai.x)
        Pi.v2pt_theory(P, A, Ai)

        # Create a new particle of mass m[i] at this point
        Pai = mechanics.Particle('Pa' + str(i), Pi, m[i])
        particles.append(Pai)

        # Set forces & compute kinematic ODE
        forces.append((Pi, m[i] * g * A.x))
        kinetic_odes.append(q[i].diff(t) - u[i])

        P = Pi

    # Generate equations of motion
    KM = mechanics.KanesMethod(A, q_ind=q, u_ind=u, kd_eqs=kinetic_odes)
    fr, fr_star = KM.kanes_equations(particles, forces)
    
    #-----------------------------------------------------
    # Step 3: numerically evaluate equations and integrate

    # initial positions and velocities – assumed to be given in degrees
    y0 = np.deg2rad(np.concatenate([np.broadcast_to(initial_positions, n), np.broadcast_to(initial_velocities, n)]))
        
    # lengths and masses
    lengths = np.broadcast_to(lengths, n)
    masses = np.broadcast_to(masses, n)

    # Fixed parameters: gravitational constant, lengths, and masses
    parameters = [g] + list(l) + list(m)
    parameter_vals = [9.81] + list(lengths) + list(masses)

    # define symbols for unknown parameters
    unknowns = [Dummy() for i in q + u]
    unknown_dict = dict(zip(q + u, unknowns))
    kds = KM.kindiffdict()

    # substitute unknown symbols for qdot terms
    mm_sym = KM.mass_matrix_full.subs(kds).subs(unknown_dict)
    fo_sym = KM.forcing_full.subs(kds).subs(unknown_dict)

    # create functions for numerical calculation 
    mm_func = lambdify(unknowns + parameters, mm_sym)
    fo_func = lambdify(unknowns + parameters, fo_sym)

    # function which computes the derivatives of parameters
    def gradient(y, t, args):
        vals = np.concatenate((y, args))
        sol = np.linalg.solve(mm_func(*vals), fo_func(*vals))
        return np.array(sol).T[0]

    # ODE integration
    return odeint(gradient, y0, times, args=(parameter_vals,))''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('get_xy_coords')

st.code('''def get_xy_coords(p, lengths):
    """
    Get x and y coordinates from lengths of pendulums and p

    Args:
        p (): points p
        lengths (list): list of lengths of each pendulum

    Returns:
        x (): numpy cumsum with all x values
        y (): numpy cumsum with all y values
    """
    # calculate x and y coordinates from lengths and points; safe values in x and y array for animation
    p = np.atleast_2d(p)
    n = p.shape[1] // 2
    lengths_neg = [ -x for x in lengths]
    zeros = np.zeros(p.shape[0])[:, None]
    x = np.hstack([zeros, lengths * np.sin(p[:, :n])])
    y = np.hstack([zeros, lengths_neg * np.cos(p[:, :n])])
    return np.cumsum(x, 1), np.cumsum(y, 1)''')

# ------------------------------------------------------------------------------------------------------------------------ #

st.subheader('animate_pendulum_multiple')

st.code('''def animate_pendulum_multiple(x, y):
    """
    Animates pendulum motion over x and y values; the figure and axes object are defined; creation of lines object;
    definition of init and animate function
    Animate function: for each frame i every line in lines is updated with values for every pendulum forthis specific index i

    Args:
        x (List): list of all x values for every pendulum
        y (list): list of all y values for every pendulum

    Returns:
        anim (): animation object
    """
    # define plot and ax; set some parameters
    title = 'Pendel motion over time (' + str(n_pendulums) + " Pendulums, " + str(n) + "-fold)"
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')
    ax.set_xlim(-1.05*sum_length, 1.05*sum_length)
    ax.set_ylim(-1.05*sum_length, 1.05*sum_length)
    ax.set_title(title)
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    ax.grid(1)

    # initiate the lines object here; lines object contains n_pendulums as line objects
    lines = [plt.plot([], [], 'o-', lw=2)[0] for _ in range(n_pendulums)]

    # initialize the lines object; iterate over every line
    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    # animate the position of each line for every frame i; iterate over every line
    def animate(i):
        for j,line in enumerate(lines):
            x = x_array[j]
            y = y_array[j]
            line.set_data(x[i], y[i])
        return lines

    # create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(times), interval=1000 * times.max() / len(times), blit=True, init_func=init)

    # closes the static plot which is not used here
    plt.close(fig)
    
    return anim''')



