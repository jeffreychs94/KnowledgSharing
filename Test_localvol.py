
import os
from Quant.Misc import *  
from Quant.LocalVolatility import *
from Quant.FiniteDifference import *
from datetime import datetime
import numpy as np

temp_EquitySpot = 48.8
temp_lstMatDate = [0.1,0.2,0.3,0.5,366/365,731/365]
temp_lstZeroRate = [0.03,0.03,0.03,0.03,0.03,0.03]
temp_arrlst_IVvalue = np.divide([[48.390,45.674,43.874,43.024,42.707],[42.699,41.131,39.959,39.399,38.841],[41.294,40.027,39.107,38.655,38.159],[38.082,37.264,36.621,36.265,35.961]],100)
temp_lst_IVdate = [30/365,91/365,182/365,366/365]
temp_arrlst_IVstrike = [[0.8,0.9,1,1.1,1.2],[0.8,0.9,1,1.1,1.2],[0.8,0.9,1,1.1,1.2]]
temp_lst_IVstrike = [0.9,0.95,1,1.05,1.10]


def load_localvol(filename = None, spot = temp_EquitySpot, lst_curveMat = temp_lstMatDate , lst_curveZero = temp_lstZeroRate, lst_IVMat = temp_lst_IVdate, arrlst_IVstrike = temp_arrlst_IVstrike, arrlst_IVvalue = temp_arrlst_IVvalue  ):

    equity = Equity.create(spot = spot, lst_curveMat = temp_lstMatDate , lst_curveZero = temp_lstZeroRate, lst_IVMat = temp_lst_IVdate, arrlst_IVstrike = temp_arrlst_IVstrike, arrlst_IVvalue = temp_arrlst_IVvalue  )
    ivol = ImpliedVol.create(equity,temp_lst_IVstrike,temp_lst_IVdate,temp_arrlst_IVvalue)
    localvol = LocalVol.from_calibration(ivol)
    
    return localvol


if __name__ == '__main__':
    np.set_printoptions(precision=14)
    localvol_equity = load_localvol()
    

