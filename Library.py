import sympy as smp
import numpy as np
import itertools as it
from GeneralRelativity.Objects import Functions

class Tensorlib:
    def __init__(self, Coordinates = [],
                       Parametrised = False, 
                       Parameter = str):
        C = Coordinates
        if Parametrised: 
            P = smp.symbols('{}'.format(Parameter))
            C[0], C[1], C[2], C[3] = smp.symbols('{} {} {} {}'.format(C[0],C[1],C[2],C[3]), cls = smp.Function)
            C[0], C[1], C[2], C[3] = C[0](P), C[1](P), C[2](P), C[3](P)
            self.Coordinates = C
            
        if not Parametrised:
            C[0], C[1], C[2], C[3] = smp.symbols('{} {} {} {}'.format(C[0],C[1],C[2],C[3]))
            self.Coordinates = C
        
    def gensym(self):
        C = self.Coordinates
        A , B = smp.symbols('A B', cls = smp.Function)
        A = A(C[1])
        B = B(C[1])
        T = smp.MutableDenseNDimArray([[-A,0,0,0],[0,B,0,0],[0,0,C[1]**2,0],[0,0,0,C[1]**2*smp.sin(C[2])**2]])
        return T
       
    def EM_strssen(self, Metric, Electromagnetic_Tensor):
        G = Metric
        F = Electromagnetic_Tensor
        Ginv = Functions(Metric, self.Coordinates, 4).Ginv()
    
    
        F_ud = smp.MutableDenseNDimArray(smp.zeros(4**2),(4,4))
        for i in range(4):
            for j in range(4):
                for d in range(4):
                    F_ud[i,j] = Ginv[i,d]*F[d,j]
        F1 = F_ud
    
        F_uu = smp.MutableDenseNDimArray(smp.zeros(4**2),(4,4))
        for i in range(4):
            for j in range(4):
                for d in range(4):
                    for s in range(4):
                        F_uu[i,j] = Ginv[i,d]*F[d,s]*Ginv[s,j]
        F2 = F_uu
     
        T = smp.MutableDenseNDimArray(smp.zeros(4**2),(4,4))
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    for s in range(4):
                        for d in range(4):
                            T[i,j] += F[i,d]*F2[d,s] - (1/4)*G[i,j]*F[k,s]*F1[k,s]
        return T

    
    def Minkowsky(Polar = True, Cartesian = False):
        if Cartesian:
           return smp.MutableDenseNDimArray([[-1,0,0,0],[0,1,0,0],[0,0,1,0], [0,0,0,1]])
        if Polar:
           t, r, theta, phi = smp.symbols('t r theta phi')
           return smp.MutableDenseNDimArray([[-1,0,0,0],[0,1,0,0],[0,0,r**2,0], [0,0,0,r**2*smp.sin(theta)**2]])
 
 
    def Electromagnetic_Tensor(E,M):
        return smp.MutableDenseNDimArray([[0,E[0],E[1],E[2]],[-M[0],0,-M[2],M[1]],[-M[1],M[2],0,-M[0]], [-E[2],-M[1],M[0],0]])
 
    
    def Reissner_Nordstrom(self,S_radius, Charge):
        C = self.Coordinates
        P = S_radius
        Q = Charge
        
        return smp.MutableDenseNDimArray([[-(1 - (P)/(C[1]) + (Q**2)/(C[1]**2)),0,0,0],[0,1/(1 - (P)/(C[1])  + (Q**2)/(C[1]**2)),0,0],[0,0,C[1]**2,0], [0,0,0,C[1]**2*smp.sin(C[2])**2]])

    def Schwarzschild(self, S_radius):
        
        C = self.Coordinates
        S_radius = smp.symbols('{}'.format(S_radius))
        P = S_radius
 
        return smp.MutableDenseNDimArray([[-(1 - (P)/(C[1])),0,0,0],[0,1/(1 - (P)/(C[1])),0,0],[0,0,C[1]**2,0],[0,0,0,C[1]**2*smp.sin(C[2])**2]])
       
    def Schwarzschild_AdS(self, S_Radius, C_Radius):
        C = self.Coordinates
        S_Radius[0] = smp.symbols('{}'.format(str(S_Radius[0])))
        C_Radius[0] = smp.symbols('{}'.format(str(C_Radius[0])))
        
        return smp.MutableDenseNDimArray([[-(1 + (C[1]**2)/(C_Radius**2) - (S_Radius)/(C[1])),0,0,0],[0,1/(1 + (C[1]**2)/(C_Radius**2) - (S_Radius)/(C[1])),0,0],[0,0,C[1]**2,0],[0,0,0,C[1]**2*smp.sin(C[2])**2]])
      
      
    def Kerr(self,S_radius, Charge):
        C = self.Coordinates
        P = S_radius
        Q = Charge
        
        return smp.MutableDenseNDimArray([[-(1 - (P)/(C[1]) + (Q**2)/(C[1]**2)),0,0,0],[0,1/(1 - (P)/(C[1])  + (Q**2)/(C[1]**2)),0,0],[0,0,C[1]**2,0], [0,0,0,C[1]**2*smp.sin(C[2])**2]])