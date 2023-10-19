import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt

class RLC():

    """
    This class initializes the ode of series RLC circuit and
    solves it. The ode is in the form of y''+2*alpha*y'+omega^2*y=0.
    """

    def alpha(r, l):
        """
        Defines alpha parameter for the ode.

        Inputs:
        -------
        r: float, resistance of the system.
        l: float, inductor of the system.

        Returns:
        --------
        alp: float, alpha parameter of the system.
        """
        alp = r/2*l
        return alp
    
    def omega_0(l, c):
        """
        Defines omega parameter for the ode.

        Inputs:
        -------
        l: float, inductor of the system.
        c: float, capacidance of the system.

        Returns:
        --------
        omega: float, omega parameter of the system.
        """

        omega = 1/np.sqrt(l*c)
        return omega
    
    def rlc_ode(t, v, r, l, c):
        """
        Defines the ode with vector v = [y', y].
        
        Inputs:
        -------
        t: list, time span for limits of integration.
        v: list, vector that defines y' and y.
        r: float, resistance of the system.
        l: float, inductor of the system.
        c: float, capacidance of the system.

        Return:
        -------
        ode: list, ode in the form of
        y'' = -2*alpha*y' - 2*omega^2*y        

        """
        alp = alpha(r, l)
        omeg = omega_0(l, c)
        ode =  [-2*alp*v[1], -2*omeg**2*v[0]]
        return ode
    
    def solve_ivp(r, l, c, v0, t):
        """
        Solves the ode with scipy.integrate.solve_ivp.

        Inputs:
        -------
        r: float, resistance of the system.
        l: float, inductor of the system.
        c: float, capacidance of the system.
        v0: list, initial conditions.   
        t: list, time span for limits of integration.

        Returns:
        --------
        y: numpy_array, solutions for y
        time: numpy array, time steps.
        """
        soln = si.solve_ivp(rlc_ode, t, v0, args=(r, l, c))
        y = soln.y
        time = soln.t
        return y, time
