import numpy as np
from scipy import special as sf

class Photon():
    
    def __init__(self, T , N):
        self.T = T #Temperature in Kelvin (5800 is approximately the BB-temp of sun)
        self.N = N #number of samples we want
        self.f = [] #list of frequencies, to be filled
        self.c= 3*10**8 #constant speed of light in m/s
        self.h = 6.6261*10**(-34) #Planck's constant in J/s
        self.kB = 1.381*10**(-23) #Boltzmann constant in 1/K
        self.wF = 3.157*self.kB*self.T/self.h #This is Wien's frequency, using Wien's law, calculate the peak frequency, use this to get the domain of interest
        
    def p_law(self, f):
        """"
        This function calculates the result of
        the Plank's Law.

        Parameters
        ----------
        f : float
            Frequency value of the photons.

        Returns
        -------
        float
            Result of the Plank's Law.

        """
        return (2*self.h*f**3/(self.c**2))*(1/(np.exp(self.h*f/(self.kB*self.T))-1))

    def random_frequency(self, N=None):
        """"
        This function generates the random frequencies. 

        Parameters
        ----------
        N : int, optional
            Number of samples. The default is None.

        Returns
        -------
        f : int
            frequency values.
        """
        f=[]
        if N is None:
            N=self.N
            
        i=0
        while i<N:
            nu = np.random.uniform(low = .001*self.wF, high = 5*self.wF)
            B = np.random.uniform(high = 1.05*self.p_law(self.wF))
            if B<self.p_law(nu):
                f.append(nu)
                i+=1
     
        return (f)
    
    def coupler(self, a, sigma, N = None):
        """
        Parameters
        ----------
        a : float
            side length of the circuit
        sigma : float
            standart deviation
        N : int, optional
            Number of samples. The default is None.

        Returns
        -------
        epc : TYPE
            DESCRIPTION.

        """
        xi=sf.erf(a/(2*sigma))
        if N is None:
            N= self.N

        u = np.random.uniform(low = 0, high = 1, size = (2,N))

        x=sigma*sf.erfinv(2*xi*u[0]-xi)   
        y=sigma*sf.erfinv(2*xi*u[1]-xi)
        
        theta = np.random.uniform(high = np.pi/2, size= N)
        
        epc=(x**2+y**2+a**2)/(2*a**2)*(np.cos(theta)-np.sin(theta))
        
        return (epc)
        
        
        
        
        
        
        