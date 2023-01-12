   
from math import log
import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve 
from Quant.Misc import *

class FDM1D:
    
    __theta__ = 0.5 # 0.5 for Crank-Nicolson
    
    """PDE: dU/dt = a * d2U/dx2 + b * dU/dx + c * U"""
    
    def __init__(self, stdev, *, fwd=1.0, nt=300, nx=300, ns= 3, fn_abc=None, bound=(-np.inf, np.inf),spot):      
        
        self.nt = nt        
        self.bound = bound
       
        xl = np.log(max(bound[0], spot*np.exp(-ns*stdev))) #spot step upper or lower #stdevis already max of BS
        xu = np.log(min(bound[1], spot*np.exp(ns*stdev)))
        
        self.xv = HashableArray(np.linspace(xl, xu, max(nx,50) + 1)) #vector for spot step #adjusted for AM vector
        self.dx = (self.xv[-1] - self.xv[0]) / (self.xv.size - 1) #spot step #dx
        #logtxt('Spot Step Size : ' + str(self.dx))
        self.fn_abc = fn_abc #pde coefficient function receive
        ###abc stands for the coefficient?
        self.eye = sparse.eye(self.xv.size) #identity matrix

    def __set_vector_and_matrix(self, t, yv,bool_vectorOutput = 0,bool_Output = 0):
        
        if bool_vectorOutput == 1 and t > 0.26 and t <0.27  : 
            input('Checkpoint 4')
            print(t)
            print(self.xv)



        a, b, c = self.fn_abc(t, self.xv,bool_vectorOutput) #coefficients used? #xv is vector for spot
        a /= self.dx ** 2 # for 2nd derivative #These are all the coeffcients
        b /= self.dx * 2

        if bool_vectorOutput == 1 and t > 0.26 and t <0.27  : 
            input('Checkpoint 5')
            print(a)
            



        # print(a)
        # input()
        dm = (a - b)[1:] #Lower Diagonal
        d0 = c - a * 2                
        dp = (a + b)[:-1]

        #Null Gamma Condition
        d0[0] = c[0] - b[0] * 2 
        dp[0] = b[0] * 2   
    
        d0[-1] = c[-1] + b[-1] * 2
        dm[-1] = -b[-1] * 2

        return yv, sparse.diags((dm,d0,dp), (-1,0,1)) #Creating triagonal matrix


    def evolve(self, start, end, yv, x=None, bool_vectorOutput = 0): # forward PDE if dt > 0 otherwise backward PDE 
        
        #for np.searchsorted
        if start > end:
            start *= 1 - 1e-14
            end *= 1 + 1e-14
        else:
            start *= 1 + 1e-14
            end *= 1 - 1e-14

        tv = np.linspace(float(start), float(end), max(abs(start - end) * self.nt, 18) + 1)

    
        dt = (tv[-1] - tv[0]) / (tv.size - 1)

        wp = dt * FDM1D.__theta__   #Crank nicholson, hence in between

        wn = dt - wp #Crank nicholson mid step

        if start > 0.082 and start <0.092 and bool_vectorOutput == 1 :
            print('This is tv, wp and dt')
            print(tv)
            print(wp)
            print(dt)     

        
        if bool_vectorOutput == 2:
            print(tv)
            print(wp)
            print(yv)
            input('CheckPoint 9')
    
        yv, mm = self.__set_vector_and_matrix(start, yv, bool_vectorOutput = 2) #yv doesn't change here
    
        if start > 0.082 and start <0.092 and bool_vectorOutput == 1 :
            input('Checkpoint 3')
            print(mm)
            print(yv)
            input('Checkpoint 3')
            print(pd.DataFrame.sparse.from_spmatrix(mm))
            input('Checkpoint 3 mm')

        if start > 0.091 and start <0.092 and bool_vectorOutput == 1 :
            print('This is start {}, and this is end {}'.format(start,end))
            print('Check Here')
            print(yv)
            
            input() 
        
        for t in tv[1:]:

            b, mm = self.__set_vector_and_matrix(t, yv + mm * yv * wp, bool_vectorOutput) #yv only changes by input. b is whatever the input there is

            yv = spsolve(self.eye - mm * wn, b)

            if bool_vectorOutput == 1 and t > 0.091 and t <0.092 :
                print('Checkpoint 2')
                print(yv)
                #input()

            
            if bool_vectorOutput == 2:
                print('Checkpoint 7')
                print(t)
                print(yv)
                input()

        return yv

    def functionize(self, yv):

        return LinearFlat(self.xv, yv)

