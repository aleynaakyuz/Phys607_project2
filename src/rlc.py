import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt

class RLC():

    """
    This class initializes the ode of series RLC circuit and
    solves it. The ode is in the form of y''+2*alpha*y'+omega^2*y=0.
    """
    
    def __init__(self, r,l,c, tint, ic):
        self.alpha = r/(2*l)
        self.omega = 1/np.sqrt(l*c)
        self.tint = tint
        self.ic = ic

    
    def odes(self, t, x):
        dx = np.zeros((2))
        dx[0] = x[1]
        dx[1] = -self.alpha*x[1]-self.omega**2*x[0]
        return dx
    
    def power(self, r, x):
        I_rms = np.sqrt(np.mean(x**2))
        return I_rms**2*r
        
    def solve_ivp(self, n):
        '''
        This function puses the scipy library to compute the 4th order Runge-Kutta method integration to find the posistion of the oscillator as a function of time

        Parameters
        ----------
        n : int
            number of points to analyze inbetween the initial and final time

        Returns
        -------
        t : ndarray
            time series used to calculate x
        x : ndarray
            posistion series calculated with t

        '''
        step = (self.tint[1]-self.tint[0])/n

        Xres = si.solve_ivp(
            self.odes, self.tint, self.ic, method='RK45', max_step=step)

        return Xres.t, Xres.y[0], Xres.y[1]
