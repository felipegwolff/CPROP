'''
  CPROP  - Converge Properties
  
  ATENTION: The Poperties are Pogressive wrong in the superheated region
  if you are working in this region please verify the Poperties in the range.
  
  Converge Engineering Solutions: www.converge-es.com
  
  Authors:
    Fabio C. Canesin <canesin@converge-es.com> 03/19/2012
    
'''
import math

#----------------------- Numerical methods needed in the library---------------
class NRB:
    '''
    Root finding using Newton-Raphson-Bisection, the method should safe and fast.
    IE: tends to quadratic convergence and to be unconditionally stable.
    
    USAGE:
    
        >>> import CPROP
        >>> NR = CPROP.NRB(func, dfunc, xmin, xmax, TOL=1e-8)
        >>> NR.solve()

    '''

    def __init__(self, func, dfunc, xmin, xmax, TOL=1e-8):
        self.func = func
        self.dfunc = dfunc
        self.a = xmin
        self.b = xmax
        self.TOL = TOL
        
    def solve(self):
        fa = self.func(self.a)
        if (fa == 0.0): return self.a
        fb = self.func(self.b)
        if (fb == 0.0): return self.b
        x = 0.5*(self.a+self.b)
        count = 0
        while count < 1001:
            fx = self.func(x)
            if abs(fx) < self.TOL: return x
            #Reduce brackets around root
            if fa*fx < 0.0:
                self.b = x
            else:
                self.a = x
            #Try a Newton-Raphson step
            dfx = self.dfunc(x)
            try: dx = - fx/float(dfx)
            except ZeroDivisionError: dx = self.b-self.a
            x = x + dx
            #If the result is outside the brackets, use a bisection step
            if (self.b - x)*(x - self.a) < 0.0:
                dx = 0.5*(self.b-self.a)
                x = self.a + dx
            #Check convergence
            if abs(dx) < self.TOL*max(abs(self.b),1.0): return x
            count +=1
        raise ValueError('Unable to solve equation :( .. sorry!')


#---------------------- EQUATION OF STATE (EOS) CLASSES ----------------------------

class MartinHou(NRB):
    '''
    Martin-Hou equation of state class. Should be used by fluid classes
         
    '''

    def __init__(self):
        '''
        Initiate MAH constants for all methods frhom selected fluid, to add a new
        fluid just copy and paste the "if" condition and replace the constants. The
        DuPond site Povides datasheets for the fluid with all the needed constants.
        '''

        ## Initiate MAH constants for all methods, see R134a for reference.

    def Cv0(self, T):
        '''
        Calculate the ideal Gas Heat Capacity (at constant volume), as a function
        of temperature.

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.Cv0(T)

        '''

        return self.cva+self.cvb*T+self.cvc*(T**2)+self.cvd*(T**3)+(self.cvf/(T**2))
      
    def Psat(self, T):
        '''
        Calculate the Vapor Pessure as a function of temperature (T [K]).

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.Psat(T)

        '''
        return 10**(self.pA + self.pB/T + self.pC*math.log10(T) + self.pD*T + \
        self.pE*((self.pF-T)/T)*math.log10(self.pF-T))

    def rhof(self, T):
        '''
        Calculate the Density of the Saturated Liquid as a function of temperature (T [K]).

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.rhof(T)

        '''
        Tr = T/self.Tc
        try:
            return  self.Af + self.Bf*(1-Tr)**(0.333333)+self.Cf*(1-Tr)**(0.666666)+ \
                    self.Df*(1-Tr) + self.Ef*(1-Tr)**(1.333333)
        except ValueError:
            raise ValueError('Temperature is above critical temperature')

    def P(self, T, v):
        '''
        Calculate the Pessure (P [kPa]) frhom the temperature (T [K]) and
        specific volume (v [m3/kg]).

        USAGE:
        
        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.P(T, v)
        
        '''
        
        z=v-self.Pb
        Tr = T/self.Tc
        return self.R*T/z + \
                (self.A[0]+self.B[0]*T+self.C[0]*math.exp(-self.Pk*Tr))/(z**2) + \
                (self.A[1]+self.B[1]*T+self.C[1]*math.exp(-self.Pk*Tr))/(z**3) + \
                (self.A[2]+self.B[2]*T+self.C[2]*math.exp(-self.Pk*Tr))/(z**4) + \
                (self.A[3]+self.B[3]*T+self.C[3]*math.exp(-self.Pk*Tr))/(z**5)
        
    def rho(self, T, P):
        '''
        Calculate the density (rho [kg/m3]) frhom the temperature (T [K])
        and Pessure (P [kPa]).

        The density is calculated iteratively using the Newton-Raphson method. 
        A ValueErrhor is raised if the method fails to converge.

        USAGE:
        
        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.rho(T, P)

        '''
        Tr = T/self.Tc

        ## z = (1/rho-b)
        ## Define f(z) = 0
        def frho(z):
            'f(z) function'
            return  -P + self.R*T/z + \
                    (self.A[0]+self.B[0]*T+self.C[0]*math.exp(-self.Pk*Tr))/(z**2) + \
                    (self.A[1]+self.B[1]*T+self.C[1]*math.exp(-self.Pk*Tr))/(z**3) + \
                    (self.A[2]+self.B[2]*T+self.C[2]*math.exp(-self.Pk*Tr))/(z**4) + \
                    (self.A[3]+self.B[3]*T+self.C[3]*math.exp(-self.Pk*Tr))/(z**5)
            
        ## Define f'(z) = 0
        def dfrho(z):
            'f(z) function'
            return -self.R*T/z**2 - \
                (2*(self.A[0]+self.B[0]*T+self.C[0]*math.exp(-self.Pk*Tr)))/z**3- \
                (3*(self.A[1]+self.B[1]*T+self.C[1]*math.exp(-self.Pk*Tr)))/z**4- \
                (4*(self.A[2]+self.B[2]*T+self.C[2]*math.exp(-self.Pk*Tr)))/z**5- \
                (5*(self.A[3]+self.B[3]*T+self.C[3]*math.exp(-self.Pk*Tr)))/z**6            

        ## Verify if it is fluid
        if(Tr <1.0 and P >= self.Psat(T)):
            return self.rhof(T)

        ## Use Newton-Raphson-Bisection method to calculate z, them recover rho
        else:
            try:
                nrb_MHrho = NRB(frho, dfrho, 1e-12, 1e12, TOL=1e-5)
                return 1/(nrb_MHrho.solve() + self.Pb)
         
            except OverflowError:
                raise ValueError('Unable to resolve with this T and P.')
            
    def s(self, T, rho):
        '''
        Calculate the enTopy (s [kJ/kg*K]) frhom the temperature (T [k])
        and density (rho [kg/m3]).

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.s(T, rho)
              
        '''
        Tr = T/self.Tc
        return  self.cvinf/(self.nexp*self.zc)*(T**self.nexp-1.)+ \
                (math.log((1./rho-self.Pb)/(1.-self.Pb)))/self.zc+(-self.b2r+ \
                self.c2r*5.475*math.exp(-self.Pk*Tr))/(1./rho-self.Pb)+(-self.b3r+ \
                self.c3r*5.475*math.exp(-5.475*T))/(2.*(1./rho-self.Pb)**2)+ \
                (-self.b5r+self.c5r*5.475*math.exp(-5.4*T))/(4.* \
                (1./rho-self.Pb)**4)

    def snorm(self, T, rho):
        '''
        Directly calculate the normalized enTopy (s-s_c)/(R*T_c)
        frhom the temperature (T) and density (rho).

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.snorm(T, rho)
                
        '''
        # Normalise with the enTopy at critical point (s_c)
        s_c = self.sc(1., 1.)       

        return self.s(T, rho) - s_c
        
#---------------------------- Transport and others properties -----------------

class Refrigerant:
    '''
    Class that defines a normal fluid transport and other properties.
    
    '''

    def sig(self, T):
        '''
        Compute the Surface Tension as function of temperature (T [K]).
        Uses the Sastri and Rao (1995) equation.

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.sig(T)
        
        '''
        SR = self.sK*((self.Pc*0.01)**self.sX)*(self.Tb**self.sY)*(self.Tc**self.sZ)
        Tr = T/self.Tc
        if (Tr > 1): raise ValueError('Sig of supercritic fluid ?')
        Tbr = self.Tb/self.Tc
        return SR*(((1-Tr)/(1-Tbr))**self.sm)

    def eta(self, T, P): #TODO!
        '''
        Model from: Model for the Viscosity and Thermal Conductivity of Refrigerants,
        Including a New Correlation for the Viscosity of R134a.
        Marcia L. Huber,* Arno Laesecke, and Richard A. Perkins
        
        Ind. Eng. Chem. Res. 2003, 42, 3163-3178

        USAGE:

        >>> from CPROP import R134a
        >>> r = R134a()
        >>> r.eta(T, P)
        
        '''

        pass

#------------------------------------ Fluids ----------------------------------

class R134a(MartinHou, Refrigerant):
    '''
    Implementation of R134 (HFC-134a) fluid using Martin-Hou equation of state.
    Inherits from Martin-Hou EoS and NormalFluid transport properties
         
    '''

    def __init__(self):
        '''
        Initiate MAH constants for all methods frhom selected fluid, to add a new
        fluid just copy and paste the "if" condition and replace the constants. The
        DuPond site Povides datasheets for the fluid with all the needed constants.
        '''
        self.kb = 1.380648813e-23 #Boltzman constant in SI [J/K]
        self.Na = 6.0221412927e-23 #Avogrado constant [mol^-1]
        ## Cnstants for R134a - From DuPont datasheet
        self.MW = 102.03 #[g/mol] - molecular weight
        self.Tb = 247.076 #[K] - Normal boiling point (at 1 atm)
        self.Tc = 374.23 #[K] - critical temperature
        self.Pc = 4060.3 #[kPa] - critical Pessure
        self.Vc =  0.00194 #[m3/kg] - critical volume
        self.rhoc = 511.9 #[kg/m3] - critical density
        self.R = 0.0815 #[kJ/kg*K] - gas constant
        self.omega = 0.32684 #Acentric factor
        self.kappa = 0.0 #Association factor
        self.sigma = 0.468932# [nm] minimum Leonard-Jones potential
        self.EpByKap = 299.363 #[K] Lennard-Jones coefficient epsilon/kappa
        self.mu = 2.058 #Dipole moment
        self.Tt = 169.85 #triple point temperature [K]
        self.hf = 200 #[kJ/kg] - reference enthalpy at 273.15K
        self.sf = 1 #[kJ/kg*K] - reference enTopy at 273.15K
        # EoS constants, for R134a it is Martin-Hou, ex: A = [A2, A3, A4, A5]
        self.A = [-8.909485e-2, -1.016882e-3, 1.778071e-5, -7.481440e-8]
        self.B = [4.408654e-5, 2.574527e-6, -4.016976e-8, 1.670285e-10]
        self.C = [-2.074834, 2.142829e-2, -2.977911e-4, 1.255922e-6]
        self.Pb = 3.755677e-4
        self.Pk = 4.599967
        #Ideal Gas Heat Capacity constants
        self.cva = 3.154856
        self.cvb = -1.656054e-2
        self.cvc = 4.353378e-5
        self.cvd = -3.754497e-8
        self.cvf = -3.023189e+4
        #Vapor Pessure constants
        self.pA = 4.069889e+1
        self.pB = -2.362540e+3
        self.pC = -1.306883e+1
        self.pD = 7.616005e-3
        self.pE = 2.342564e-1
        self.pF = 3.761111e+2
        #Density of Saturated Liquid constants
        self.Af = 5.281464e+2
        self.Bf = 7.551834e+2
        self.Cf = 1.028676e+3
        self.Df = -9.491172e+2
        self.Ef = 5.935660e+2
        #Surface Tension constants - From CGC optimization
        self.sK = 0.15835
        self.sX = 0.53
        self.sY = -1.45
        self.sZ = 1.79
        self.sm = 1.26

    # Need to improve to acound pressure effects - Implement paper in class Refrigerants
    def etaRef(self, T, P):
        Tr = T/self.Tc
        ## Verify if it is fluid
        if(Tr <1.0 and P >= self.Psat(T)):
            if((T>=216.15) and (T<=366.15)):
                return -0.0002191*((T-273.15)**3)+0.039304*((T-273.15)**2)- \
                        3.6494*(T-273.15)+267.67
        
        else:
            if((T>=290.15) and (T<=422.15)):
                return 11.021+0.038599*(T-273.15)

    # Just overide the eta from Refrigerants for R134a - Sorry guys, deadline comming :'(
    def eta(self, T, P): 
        return self.etaRef(T, P)
