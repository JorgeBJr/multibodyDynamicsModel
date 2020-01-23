from scipy.integrate import odeint
import numpy as np
from numba import jit

## The MATLAB version is "myODE_5.m"
@jit(cache = True)
def myODE_5(Q, t_spray, F, alpha, tau0, tau_w, L1, L2, L3, L_petiole, ahead, 
            abutt, bhead, bbutt, K, c, rho_head, rho_butt, rhoA, muA, g, m1, 
            m2, echead, ecbutt, I1, I2, S_head, S_butt, betaR, tsExp):
    #Unpack the state vector (THE ORDER MATTERS)
    x, y, theta, phi, xd, yd, thetad, phid = Q 
    
    #Reynolds number calculation:
    Re_head = rhoA*(np.sqrt((xd**2)+(yd**2)))*(2*bhead)/muA #dimensionless number
    Re_butt = rhoA*(np.sqrt((xd**2)+(yd**2)))*(2*bbutt)/muA #dimensionless number
                
    #Coefficient of drag stuff:
    Cd_head = 24/np.abs(Re_head) + 6/(1 + np.sqrt(np.abs(Re_head))) + 0.4
    Cd_butt = 24/np.abs(Re_butt) + 6/(1 + np.sqrt(np.abs(Re_butt))) + 0.4
    
    #These are the coefficients for our acceleration equations imported from Mathematica 
        #(and MATLAB). Careful, this'll get messy.
    h1 = m1 + m2
    h2 = (-1)*L1*m1*np.sin(theta)
    h3 = (-1)*L2*m2*np.sin(phi)
    h4 = L1*m1*np.cos(theta)
    h5 = L2*m2*np.cos(phi)
    h6 = ((-1)*F*np.cos(alpha+theta)+(1/2)*Cd_butt*rhoA*S_butt*np.abs(xd)*xd
          +(1/2)*Cd_head*rhoA*S_head*np.abs(xd)*xd+(-1)*L1*m1*np.cos(theta)*thetad**2
          +(-1)*L2*m2*np.cos(phi)*phid**2)
    h7 = (g*(m1+m2)+(1/2)*Cd_butt*rhoA*S_butt*np.abs(yd)*yd+(1/2)*Cd_head*rhoA
          *S_head*np.abs(yd)*yd+(-1)*L1*m1*thetad**2*np.sin(theta)
          +(-1)*F*np.sin(alpha+theta)+(-1)*L2*m2*phid**2*np.sin(phi))
    h8 = ((-1)*tau0+g*L1*m1*np.cos(theta)+(-1)*K*(((-1)*betaR+(-1)*np.pi+
           (-1)*theta+phi)**tsExp)+(-1)*c*((-1)*thetad+phid)+(-1)*F*L3*np.sin(alpha)-tau_w)
    h9 = (tau0+g*L2*m2*np.cos(phi)+K*(((-1)*betaR+(-1)*np.pi+(-1)*theta+phi)**tsExp)
    +c*((-1)*thetad+phid))
    
    h10 = I1+(L1**2)*m1
    h11 = I2+(L2**2)*m2

    #x double dot
    xdd = ((-1)*(h10*h11*h1**2+(-1)*h11*h1*h2**2
            +(-1)*h10*h1*h3**2+(-1)*h11*h1*h4**2+h3**2*h4**2
            +(-2)*h2*h3*h4*h5+(-1)*h10*h1*h5**2+h2**2*h5**2)**(-1)*(h10*h11*h1*h6
             +(-1)*h11*h4**2*h6+(-1)*h10*h5**2*h6+h11*h2*h4*h7+h10*h3*h5*h7
             +(-1)*h11*h1*h2*h8+(-1)*h3*h4*h5*h8+h2*h5**2*h8
             +(-1)*h10*h1*h3*h9+h3*h4**2*h9+(-1)*h2*h4*h5*h9))
  
    #y double dot
    ydd = ((-1)*((-1)*h10*h11*h1**2+h11*h1*h2**2+h10*h1*h3**2
            +h11*h1*h4**2+(-1)*h3**2*h4**2+2*h2*h3*h4*h5
            +h10*h1*h5**2+(-1)*h2**2*h5**2)**(-1)*((-1)*h11*h2*h4*h6
                          +(-1)*h10*h3*h5*h6+(-1)*h10*h11*h1*h7+h11*h2**2*h7
                          +h10*h3**2*h7+h11*h1*h4*h8+(-1)*h3**2*h4*h8
                          +h2*h3*h5*h8+h2*h3*h4*h9+h10*h1*h5*h9+(-1)*h2**2*h5*h9))

    #theta double dot
    thetadd = ((-1)*((-1)*h10*h11*h1**2+h11*h1*h2**2+h10*h1*h3**2
                +h11*h1*h4**2+(-1)*h3**2*h4**2+2*h2*h3*h4*h5
                +h10*h1*h5**2+(-1)*h2**2*h5**2)**(-1)*(h11*h1*h2*h6
                              +h3*h4*h5*h6+(-1)*h2*h5**2*h6+h11*h1*h4*h7
                              +(-1)*h3**2*h4*h7+h2*h3*h5*h7+(-1)*h11*h1**2*h8
                              +h1*h3**2*h8+h1*h5**2*h8+(-1)*h1*h2*h3*h9
                              +(-1)*h1*h4*h5*h9))

    #phi double dot
    phidd = ((-1)*((-1)*h10*h11*h1**2+h11*h1*h2**2+h10*h1*h3**2
              +h11*h1*h4**2+(-1)*h3**2*h4**2+2*h2*h3*h4*h5
              +h10*h1*h5**2+(-1)*h2**2*h5**2)**(-1)*(h10*h1*h3*h6
                            +(-1)*h3*h4**2*h6+h2*h4*h5*h6+h2*h3*h4*h7
                            +h10*h1*h5*h7+(-1)*h2**2*h5*h7+(-1)*h1*h2*h3*h8
                            +(-1)*h1*h4*h5*h8+(-1)*h10*h1**2*h9+h1*h2**2*h9
                            +h1*h4**2*h9))

    return(np.array([xd, yd, thetad, phid, xdd, ydd, thetadd, phidd]))

    #Generate simulations (in parallelized cores)
def generateSimulations(simNum, t_spray, q0, F_array, alpha_array, tau0_array, 
                        tau_w_array, timestep, RecedeFrac, costCoeff, 
                        goalCriteria, globalList):
    F = F_array[simNum]
    alpha = alpha_array[simNum]
    tau0 = tau0_array[simNum]
    tau_w = tau_w_array[simNum]
            
    Q = odeint(myODE_5, q0, t_spray, args = (F, alpha, tau0, tau_w, *globalList))
    
    #Cost function value
    cost = (costCoeff[3]*((Q[-1,4]-goalCriteria[3])**2)
    +costCoeff[4]*((Q[-1,5]-goalCriteria[4])**2)
    +costCoeff[5]*((Q[-1,6]-goalCriteria[5])**2)
    +costCoeff[0]*((Q[-1,0]-goalCriteria[0])**2)
    +costCoeff[1]*((Q[-1,1]-goalCriteria[1])**2)
    +costCoeff[2]*((Q[-1,2]-goalCriteria[2])**2))

    #Calculate the new initial conditions
    NewIC_index = int(timestep*RecedeFrac)
    NewICs = np.array([Q[NewIC_index,0], Q[NewIC_index,1], Q[NewIC_index,2], 
                       Q[NewIC_index,3], Q[NewIC_index,4], Q[NewIC_index,5], 
                       Q[NewIC_index,6], Q[NewIC_index,7]])
    
    #This will check if the applied efforts are saved correctly
    check = np.array([F, alpha, tau0, tau_w])
    
    #The list which contains our generated simulation array, the cost function 
        #value, the new initial conditions, and our check.
    PartPath = [Q, cost, NewICs, check]
   
    return(PartPath)
    #End of numOfTrajectories loop (j) (for solving ODEs)