# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 15:41:44 2016

@author: ramandeepbedi
"""
from matplotlib import pyplot as plt
import numpy as np

class SinglePendulum():
    
    def __init__(self, D, theta_0, tMin, tMax, h=0.01):
        self.D = D
        self.theta_0 = theta_0
        self.h = float(h)
        self.t = np.arange(tMin, tMax+1, h)
        if(D<0 or h<=0 or tMin>tMax):
            raise Exception()
       
        k = np.sqrt(1-(self.D**2/4))
        self.trueTheta = np.exp(-self.D*self.t/2) * ((self.theta_0 * np.cos(k*self.t)) + (self.D*self.theta_0/(2*k))*np.sin(k*self.t))
        trueDerivative = (-(self.D/2)*self.trueTheta) + np.exp(-self.D*self.t/2)*(-k*self.theta_0*np.sin(k*self.t) + (self.D*self.theta_0/(2*k))*np.cos(k*self.t))
        self.trueEnergy = (trueDerivative**2)/2 + (self.trueTheta**2)/2  
    
    def trueSolution(self, theta):
        t = self.t
        k = np.sqrt(1-(self.D**2/4))
        
        trueU = np.exp(-self.D*t/2) * ((theta * np.cos(k*t)) + (self.D*theta/(2*k))*np.sin(k*t))
        trueDerivative = (-(self.D/2)*trueU) + np.exp(-self.D*t/2)*(-k*theta*np.sin(k*t) + (self.D*theta/(2*k))*np.cos(k*t))
        trueE = (trueDerivative**2)/2 + (trueU**2)/2   
      
        
        return [trueU, trueE]
        
    def stabilityAnalysis(self, v, u):
        
        energy = []
        for i in range(len(self.t)):
            energy.append(((v[i]**2)+(u[i]**2))/2)
            
        energyChange = [0] #Required to be 0 so that the number of x coords match the number of y coords
        for i in range(1, len(self.t)):
            energyChange.append(energy[i]-energy[i-1])
        
        plt.figure()
        plt.plot(self.t, energy, 'r')   
        plt.title('Energy-Time')
        plt.xlabel('Time')
        plt.ylabel('Energy')
        plt.grid() 
        
        plt.figure()
        plt.plot(self.t, energyChange, 'y')
        plt.title('Energy Change')
        plt.xlabel('Time')
        plt.ylabel('Energy')
        plt.grid() 
        
    
    def explicitEuler(self):
     
       # v = np.ones(len(self.t))
        v = np.full(len(self.t), 0.01)
        theta = np.full(len(self.t), 0.01)
        
        for i in range(1, len(self.t)):
            v[i] = v[i-1] - (theta[i-1]+(self.D*v[i-1]))*self.h
            theta[i] = theta[i-1] + v[i-1]*self.h
            
        plt.figure()
        plt.title("Explicit Euler")
        plt.xlabel("Time / s")
        plt.ylabel("Theta - θ / °)")
        plt.plot(self.t, theta, label="Explicit Euler with h="+str(self.h))
        plt.plot(self.t, self.trueTheta, '--', label="True Solution")
        plt.grid()       
        plt.legend()
        plt.show()
        
        self.stabilityAnalysis(v, theta)
    
    def leapfrog(self):
        
        v = np.full(len(self.t), 0.01)
        theta = np.full(len(self.t), 0.01)
    
        #Euler forward used to estimate our initial (n-1)th term.    
        v[1] = v[0] - (theta[0]+self.D*v[0])*self.h    
        theta[1] = theta[0]+v[0]*self.h
        
        for i in range(2, len(self.t)):
            v[i] = v[i-2] + (2*-(theta[i-1]+self.D*v[i-1]))*self.h
            theta[i] = theta[i-2] + 2*(v[i-1])*self.h
            
        plt.figure()    
        plt.grid()
        plt.title("Leapfrog")
        plt.xlabel("Time / s")
        plt.ylabel("Theta - θ / °")
        plt.plot(self.t, theta, label="Leapfrog with h="+str(self.h))
        plt.plot(self.t, self.trueTheta, '--', label="True Solution")
        plt.legend()
        plt.show()
        
        self.stabilityAnalysis(v, theta)
        
        
    def RK4(self):
        
        v = np.full(len(self.t), 0.01)
        u = np.full(len(self.t), 0.01)
       
        for i in range(1, len(self.t)):
            
            kv1 = -(self.D*v[i-1]+u[i-1])
            ku1 = v[i-1]
            v1 = v[i-1]+kv1*(self.h/2)
            u1 = u[i-1]+ku1*(self.h/2)
            
            kv2 = -(self.D*v1 + u1)
            ku2 = v1
            v2 = v[i-1]+kv2*self.h/2
            u2 = u[i-1]+ku2*self.h/2
            
            kv3 = -(self.D*v2 + u2)
            ku3 = v2
            v3 = v[i-1]+kv3*(self.h)
            u3 = u[i-1]+ku3*(self.h)
            
            kv4 = -(self.D*v3 + u3)
            ku4 = v3
            
            v[i] = v[i-1] + (self.h/6)*(kv1 + 2*kv2 + 2*kv3 + kv4)
            u[i] = u[i-1] + (self.h/6)*(ku1+2*ku2+2*ku3+ku4)
         
        plt.title("RK4")
        plt.xlabel("Time / s")
        plt.ylabel("Theta θ / °")
        plt.grid()
        plt.plot(self.t, u, label="RK4 with h="+str(self.h)+", D="+str(self.D))
        plt.plot(self.t, self.trueTheta, '--', label="True Solution")
        plt.legend()
        plt.show()
        
        self.stabilityAnalysis(v, u)
        

    def implicitEuler(self):      
        v = np.full(len(self.t), 0.01)
        u = np.full(len(self.t), 0.01)
        
        denom = 1 / (self.h**2 + self.D*self.h +1)
        
        for i in range(1, len(self.t)):
            v[i] = (-(self.h*u[i-1])+v[i-1])*denom 
            u[i] = u[i-1] + v[i]*self.h
                      
        plt.figure()
        plt.plot(self.t, u, label="Implicit Euler with h="+str(self.h))
        plt.plot(self.t, self.trueTheta, '--', label="True Solution")
        plt.title('Implicit Euler')
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.grid()
        plt.legend()
        plt.show()
        
        self.stabilityAnalysis(v, u)
 
     #initialTheta in degrees
    def AAOExplicitEuler(self, initialTheta):              
         v = np.full(len(self.t), 0.01)
         u = np.full(len(self.t),0.01)
         u[0] = initialTheta*np.pi/180
         
         for i in range(1, len(self.t)):
             v[i] = v[i-1] + (-np.sin(u[i-1]) - self.D*v[i-1])*self.h
             u[i] = u[i-1] + v[i-1]*self.h
         
         
         plt.title("Arbitary Angle Oscillation")
         plt.plot(self.t, u, label="Explicit Euler, h="+str(self.h)+", theta="+str(np.pi*initialTheta/180)[:4])
         plt.plot(self.t, self.trueSolution(initialTheta*np.pi/180)[0], '--', label="True solution for theta = "+str(np.pi*initialTheta/180)[:4])
         plt.xlabel("Time")
         plt.ylabel("Angle")
         plt.legend()
         plt.grid()
         plt.show()
         self.stabilityAnalysis(v, u)
       
        
class DoublePendulum():
    
    def __init__(self, h, R, G, maxTime):
        self.h = float(h)
        self.R = R
        self.G = G
        self.t = np.arange(0, maxTime, h)
        
    def DoubleRK4(self):
        u, p, w, v, TotalE, Ek, Epe = np.zeros(len(self.t)), np.zeros(len(self.t)), np.zeros(len(self.t)), np.zeros(len(self.t)), np.zeros(len(self.t)), np.zeros(len(self.t)), np.zeros(len(self.t))
        #Initial Coniditon
        u[0] = 0.1
        
        for i in range(1, len(self.t)):
            kv1 = (self.R+1)*u[i-1]-(self.R+1)*p[i-1]+self.G*(1-(self.R**-1))*w[i-1]-(self.G*self.R**-1)*v[i-1]
            kw1 = -(self.R+1)*u[i-1]+self.R*p[i-1]-self.G*w[i-1]
            kp1 = v[i-1]
            ku1 = w[i-1]
            
            v1 = v[i-1]+kv1*float(self.h/2)
            w1 = w[i-1]+kw1*float(self.h/2)
            p1 = p[i-1]+kp1*float(self.h/2)
            u1 = u[i-1]+ku1*float(self.h/2)
            
            kv2 = (self.R+1)*u1-(self.R+1)*p1+self.G*(1-(self.R**-1))*w1-(self.G*self.R**-1)*v1
            kw2 = -(self.R+1)*u1+self.R*p1-self.G*w1
            kp2 = v1
            ku2 = w1
            
            v2 = v[i-1]+kv2*self.h/2
            w2 = w[i-1]+kw2*self.h/2
            p2 = p[i-1]+kp2*self.h/2
            u2 = u[i-1]+ku2*self.h/2
            
            kv3 = (self.R+1)*u2-(self.R+1)*p2+self.G*(1-(self.R**-1))*w2-(self.G*self.R**-1)*v2
            kw3 = -(self.R+1)*u2+self.R*p2-self.G*w2
            kp3 = v2
            ku3 = w2
            
            v3 = v[i-1]+kv3*self.h
            w3 = w[i-1]+kw3*self.h
            p3 = p[i-1]+kp3*self.h
            u3 = u[i-1]+ku3*self.h
            
            kv4 = (self.R+1)*u3-(self.R+1)*p3+self.G*(1-(self.R**-1))*w3 -(self.G*self.R**-1)*v3
            kw4 = -(self.R+1)*u3+self.R*p3-self.G * w3
            kp4 = v3
            ku4 = w3
            
            v[i] = v[i-1] + self.h/6. * (kv1 + 2*kv2 + 2*kv3 + kv4)
            w[i] = w[i-1] + self.h/6. * (kw1 + 2*kw2 + 2*kw3 + kw4)
            p[i] = p[i-1] + self.h/6. * (kp1 + 2*kp2 + 2*kp3 + kp4)
            u[i] = u[i-1] + self.h/6. * (ku1 + 2*ku2 + 2*ku3 + ku4)
    
        for i in range(len(self.t)):
            
            PE = (u[i]**2)/2. + self.R*(u[i]**2)/2. + self.R*(p[i]**2)/2.
            KE = (1/2. * w[i]**2) + (1/2. * self.R * (2*w[i]*v[i] + w[i]**2 + v[i]**2))
            
            TotalE[i] = float(PE) + float(KE)
                
            Ek[i] = KE
            Epe[i] = PE

        plt.title("Energy - Time")
        plt.plot(self.t, TotalE, 'r', label="Total energy, R="+str(self.R)+", G="+str(self.G))
        plt.plot(self.t, Ek, 'g', label="Kinetic Energy")
        plt.plot(self.t, Epe, 'b', label="Potential Energy")
        plt.xlabel("Time")
        plt.ylabel("Scaled Energy")
        plt.ylim(0, max(TotalE)+max(TotalE)/5)
        plt.grid()
        plt.legend()
        plt.show()
        
        plt.figure()
        plt.title("RK4 - Double Pendulum - R="+str(self.R)+", G="+str(self.G)+", h="+str(self.h))
        plt.plot(self.t, u, 'b', label="Theta")
        plt.plot(self.t, p, 'r', label="Phi")
        plt.xlabel("Time")
        plt.ylabel("Angle")
        plt.grid()
        plt.legend()
        plt.show()
        
    
    
    
    
#D=0, theta0 = 0.01, tMin = 1, tMax = 20, h=0.1
'''    
pend = SinglePendulum(0.2, 0.01, 0, 100, 0.05)
pend.explicitEuler()
pend.RK4()
pend.leapfrog()
pend.implicitEuler()
pend.AAOExplicitEuler(135)

pendula = DoublePendulum(0.06, 1, 1, 100)
pendula.DoubleRK4()

pend = SinglePendulum(0, 0.01, 0, 50, 0.2)
pend.RK4()
'''

double = DoublePendulum(0.1, 0.01, 0, 100)
double.DoubleRK4()
