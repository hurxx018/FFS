# import numpy as np
from numpy import ndarray, sqrt, exp
# from binfunctions.BinningFuncWrapper import BinningFuncWrapper

"""Module for models for Single color FCS analysis
Function Names:
    A single diffusing species in 2DG PSF  >> acf2dg
    A single diffusing species in 3DG PSF  >> acf3dg

"""

def _checkx(x):
    if isinstance(x, (int, float)):
        return ndarray(x)
    elif isinstance(x, (list, tuple)):
        return ndarray(x)
    elif isinstance(x, ndarray):
        return x
    else:
        raise TypeError("The type of x is not correct.")

def acf2dg(x, paras, tsampling=1):
    """Autocorrelation function of a single diffusing species in 2DG PSF

    paras : [G0, td, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 3:
        raise ValueError("The number of elements in paras should be 3")
    return paras[0]/(1. + t/paras[1]) + paras[2]

def acf3dg(x, paras, tsampling=1):
    """Autocorrelation function of a single diffusing species in 3DG PSF

    paras : [G0, td, r^2, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 4:
        raise ValueError("The number of elements in paras should be 4")
    return paras[0]/(1. + t/paras[1])/sqrt(1. + t/paras[1]/paras[2]) + paras[3]

def acf2dgEC(x, paras, tsampling=1):
    """Autocorrelation function of a single diffusing species in 2DG PSF in the
    presence of exponential correlation (Gaussian Noise)

    paras : [G0, td, A0, T0, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 5:
        raise ValueError("The number of elements in paras should be 5")
    return paras[0]/(1. + t/paras[1]) + paras[2]*exp(-t/paras[3])+ paras[4]



def main():
    import matplotlib.pyplot as plt
    from numpy import arange
    x = (arange(100000) + 1)
    plt.plot(x*0.00001, acf2dg(x, [0.01, 0.001, 0.], 0.00001), label='2DG')
    plt.plot(x*0.00001, acf3dg(x, [0.01, 0.001, 25., 0.], 0.00001), label='3DG')
    plt.plot(x*0.00001, acf2dgEC(x, [0.01, 0.001, 0.0001, 1., 0.], 0.00001), label='2DGEC')
    plt.xscale("log")
    plt.legend()
    plt.show()

if __name__=="__main__":
    main()
