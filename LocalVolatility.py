from math import log, exp, sqrt
from telnetlib import TN3270E
import numpy as np
from scipy.optimize import minimize 
from scipy.stats import norm 
from Quant.Misc import * 
from Quant.FiniteDifference import *
import pandas 

def ascending_pairs(*args):

    steps = np.sort(np.unique(np.hstack(args)))

    return zip(steps[:-1], steps[1:])

class Curve:

    def __init__(self,t0, fn):
        self.t0 = t0
        self.__zeroRate= fn

    def __call__(self, t): # t can be float or numpy.arrays

        return self.__zeroRate(t)

    @classmethod
    def create(cls,lst_MatDate,lst_Zero):

        t0 = 0

        tenors = np.array(lst_MatDate)

        zero = np.array(lst_Zero)

        return cls(t0,LinearFlat(tenors,zero))


class Equity:

    def __init__(self, tradedate, spot, zero, div = None):
        self.t0 = tradedate
        self.spot = spot
        self.zero = zero

    @classmethod
    def create(cls, spot, lst_curveMat, lst_curveZero, lst_IVMat, arrlst_IVstrike, arrlst_IVvalue):
        t0 = 0
        zero = Curve.create(lst_curveMat, lst_curveZero)
        return cls(t0,spot,zero)
   
    def forward(self, t): 
        return self.spot / exp(-self.zero(t)*t)


class ImpliedVol:
    class Smile(LinearExtra):          
        def __init__(self, strikes, vols):
            super().__init__(strikes, vols)

        def __call__(self,k): 
            return super().__call__(k)

    def __init__(self, equity):
        self.equity = equity                         
        self.tenors = None
        self.smiles = None
        self.maxvol = None
            
    def set(self, t, strikes, vol):
        """update implied volatility smile of the longest maturity""" 
        smile = self.Smile(strikes, vol)
        if self.tenors is None: # initiating
            self.tenors = np.array([t])
            self.smiles = np.array([smile])
            self.maxvol = np.array([max(vol)])
        elif t > self.tenors[-1]: # bootstrapping, insert one more entry
            self.tenors = np.append(self.tenors, t)  
            self.smiles = np.append(self.smiles, smile)  
            self.maxvol = np.append(self.maxvol, np.array(max(vol)))
        else:
            raise BaseException('tenor must be >= current tenor structure ...')
        return self

    @classmethod
    def create(cls,equity,strikes,lst_tenors,arrlst_vol): #strike is in moneyness
        ivol = cls(equity)
        strikes = np.array(strikes) * equity.spot
        
        for i in range(0,len(lst_tenors)):
            ivol.set(lst_tenors[i], strikes, arrlst_vol[i])
        
        return ivol
        
    def __call__(self, t, z): # interpolate implied vol; not (yet) interpolated in time dimension        
        
        i = np.searchsorted(self.tenors[:-1], t)

        return self.smiles[i](z)

    def stdev(self, t, k = None): # standard deviation given expiry and log-moneyness 
        
        if k is None: #return ATM
            k = self.equity.spot

        return self(t, k) * (t) ** 0.5  

    def getmaxvol(self,t): ###Temporary Solution, Need Upgrade for >1 Pillar
        i = np.searchsorted(self.tenors[:-1], t)
        return self.maxvol[i]

    def expiries(self, start, end): # tenors between start and end time
        return self.tenors[(start < self.tenors) & (self.tenors < end)]

    def norm_call(self,t,k,fwd):
        
        stdev = self.stdev(t, k)
        d1 = (np.log(fwd/k)+0.5*stdev**2)/stdev
        d2 = d1 - stdev
        df = exp(-self.equity.zero(t)*t)

        return (fwd*norm.cdf(d1) - k*norm.cdf(d2))*df #Currently no discount factor


class LocalVol:

    class Smile(LinearExtra):

        grid = np.linspace(-1.5,1.5,7) # 7 strike pillers

        def __init__(self, stdev, vol, fwd):

                self.k = np.exp(stdev * self.grid) * fwd

                self.v = vol 

                super().__init__(np.log(self.k), self.v)

    class ExtendedSmile(LinearFlat):

        

        grid = np.linspace(-1.5,1.5,7) # 7 strike pillers

        def __init__(self, stdev, vol, fwd):
    
            extendedGrid= np.append(self.grid,3)
            extendedGrid = np.append(-3,extendedGrid)
            self.k = np.exp(stdev * extendedGrid) * fwd

            top_vol = LinearExtra(np.log(self.k[1:3]), vol[0:2])(np.log(self.k[0]))
            
            bot_vol = LinearExtra(np.log(self.k[-3:-1]), vol[-2:])(np.log(self.k[-1]))
            self.v = np.hstack((top_vol, vol, bot_vol))
            
            print('Extended Smile Create with topvol {}, and bottom vol {}, with strike scale {} and vol points {}'.format(top_vol, bot_vol, self.k, self.v))

            super().__init__(np.log(self.k), self.v)

#################For Output#########################

    class Output:
        def __init__(self):
            self.tenors = None
            self.pillar = None
            self.vol = None
        
        def fillLV(self,t,k,v): #k is in log moneyness
            if self.tenors is None: #Initiating
                self.tenors = np.array([t], dtype=object)
                self.pillar = np.array([[k]],dtype=object)
                self.vol = np.array([v],dtype=object)
            elif t > self.tenors[-1]: # appending
                self.tenors = np.append(self.tenors, t) 
                self.pillar = np.append(self.pillar, k)
                self.vol = np.append(self.vol, v)
            else:
                raise BaseException('tenor must be >= current tenor structure ...')

#################For Output#########################

    def __init__(self, ivol):
        self.ivol = ivol                        
        self.tenors = None
        self.smiles = None
        self.outputDic = None
    
    def set(self, t, vol, fwd, bool_final = 0):
        """update volatility smile of the longest maturity""" 
        if bool_final == 1:
            smile = self.ExtendedSmile(self.ivol.stdev(t), vol, fwd)
        else:
            smile = self.Smile(self.ivol.stdev(t), vol, fwd)
        
        if self.tenors is None: # initiating
            self.tenors = np.array([t])
            self.smiles = np.array([smile])
        elif t == self.tenors[-1]: # updating
            self.smiles[-1] = smile
        elif t > self.tenors[-1]: # appending
            self.tenors = np.append(self.tenors, t) 
            self.smiles = np.append(self.smiles, smile)
        else:
            raise BaseException('tenor must be >= current tenor structure ...')      
        return self


    def __fetch__(self, i, z):
        # z fed into here is in ln k
        return self.smiles[i](z)

    def __call__(self, t, z):
        i = np.searchsorted(self.tenors[:-1], t)
        return self.__fetch__(i, z) 
    
    
    @classmethod

    def from_calibration(cls, ivol): # calibrate from implied vol surface 
        
        localvol = cls(ivol) #Initiate class Ivol
        
        t0 = localvol.ivol.equity.t0
        
        n = 7 # number of strike pillers

        soln = None

        Output_dic = localvol.Output()

        for t in ivol.tenors: #Ascending pillars
            
            ks =  localvol.ivol.equity.forward(t) * np.exp(ivol.stdev(t, localvol.ivol.equity.spot) * cls.Smile.grid)
  
            calls_iv = ivol.norm_call(t, ks, localvol.ivol.equity.forward(t))

            def find_localvol(lv):
                
                calls_lv = localvol.set(t, lv, localvol.ivol.equity.forward(t)).__pde_fwd_call(t0, t, soln)(np.log(ks)) ###v,t is the 2 boundary conditions (2 pmaturity pillars                

                return sqrt(np.mean((calls_lv - calls_iv) ** 2)) #Difference between local vol and implied vol

            lv = minimize(find_localvol, x0=ivol(t,ks), bounds=[(0,1)] * n, method='L-BFGS-B').x #options = {'gtol': 1e-9,'ftol': 1e-7} 
            
            soln = localvol.set(t, lv, localvol.ivol.equity.forward(t),bool_final = 0).__pde_fwd_call(t0, t, soln,bool_vectorOutput = 0) # update local vol

            print(soln)

            #Extend to 9 Pillars
            localvol.set(t, lv, localvol.ivol.equity.forward(t),bool_final = 1)

            Output_dic.fillLV(t,ks,lv)

            print('%s:' % round(t,4), lv)
        
        localvol.outputDic = Output_dic

        return localvol 

    def __pde_fwd_call(self, start, end, fn=None, hide_div=False, bool_vectorOutput = 0): # normalized forward (i.e. undiscounted) call                      
        
        def fn_abc(t, k,bool_vectorOutput = 0): # k = ln(K/F) # k is already in ln strike

            

            lv = self(t, k)

            if bool_vectorOutput == 1 and t > 0.26 and t <0.27  : 
                input('Checkpoint 6')
                print(t)
                print(k)
                print(lv)

            r = self.ivol.equity.zero(t) #None Forward Rate For Now
            div = 0 #Currrently Default

            a = lv ** 2 / 2 #This is the coefficient of a?           
            b = -((r-div)+a)
            c = - div

            return a, b, np.repeat(c,np.shape(a)[0]) #create vector with same size

        maxstdev = self.ivol.getmaxvol(end) *(end**0.5) ##max of stdev
        
        fdm = FDM1D(maxstdev, nt=18, nx=3, ns=3, fn_abc = fn_abc,spot = self.ivol.equity.spot)
        
        yv = np.maximum(self.ivol.equity.spot - np.exp(fdm.xv), 0) #call price vector
        
        expiries = self.ivol.expiries(start, end) #2 pillars from ascending pairs

        for p, q in ascending_pairs(start, expiries, end): 
            
            if bool_vectorOutput == 1 and end > 0.25:
                
                print('Boolean to 2 start,end {} {}'.format(start,end))
                bool_vectorOutput = 2

            yv = fdm.evolve(p, q, yv, bool_vectorOutput = bool_vectorOutput)
        
        if bool_vectorOutput == 2 :

            return yv

        else:

            return fdm.functionize(yv)
 

if __name__ == '__main__': 
    pass

