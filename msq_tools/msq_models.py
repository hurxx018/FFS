import numpy as np
from binfunctions.BinningFuncWrapper import BinningFuncWrapper

"""Module for models for MSQ analysis
Function Names:
    A single diffusing species    >> biasQestimator
    The effect of shot noise      >> biasQestimator_shotnoise
    A single diffusing species and an exponential correlation
    (due to the Gaussian Noise)   >> biasQestimator_expcor
    Two diffusing species         >> biasQestimator_Twospecies
    Three diffusing species       >> biasQestimator_threespecies
    Two diffusing species and an exponential correlation
    (due to the Gaussian Noise)   >> biasQestimator_Twospecies_expcor


"""


def biasQestimator(x, paras, qshift=1, tsampling=None):
    """Qestimator bias model for a system consisting of 3DG PSF and a single
    diffusing species
    Inputs: x, paras
    x : a numpy array of segment times
    paras : [Q0 = gamma2*lambda*Ts, td -> diffusion time, alpha -> z0^2/w0^2],
            where Ts is a sampling time

    qshift (integer) : 1 -> Q1, 0 -> Q0
    tsampling (float) : sampling time


    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if not isinstance(qshift, int):
        raise TypeError("qshift is integer")
    else:
        v = qshift
    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling
    bin2 = BinningFuncWrapper(order=2).calculate

    Q0 = paras[0]
    td = paras[1]
    alpha = paras[2]

    result = np.zeros(tx.size, dtype=np.float64)

    if v == 0:
        result += Q0*bin2(tx, alpha, td)/np.power(t, 2)
        result += -1./tx
        result += -Q0 * bin2(tx*t, alpha, td)/np.power(tx*t, 2)
    else:
        result += Q0*(0.5*bin2((v + 1)*t, alpha, td)
                        + 0.5*bin2((v - 1)*t, alpha, td)
                        - bin2(v*t, alpha, td))/np.power(t, 2)
        index = np.where(tx/2 >= v)[0]
        if index.size == tx.size:
            result += -(tx - 2*v)/np.power(tx - v, 2)
            result += -Q0 * (bin2(tx*t, alpha, td) + bin2((tx - 2*v)*t, alpha, td)
                            - 2.*bin2(v*t, alpha, td))/2./np.power((tx-v)*t, 2)
        else:
            result[index] += -(tx - 2*v)/np.power(tx[index] - v, 2)
            result[index] += -Q0 * (bin2(tx[index]*t, alpha, td) + bin2((tx[index] - 2*v)*t, alpha, td)
                            - 2.*bin2(v*t, alpha, td))/2./np.power((tx[index] - v)*t, 2)
            index1 = np.where(tx/2 < v)[0]
            result[index1] += -Q0 * (bin2(tx[index1]*t, alpha, td) + bin2((2*v + tx[index1])*t, alpha, td)
                            - 2.*bin2(v*t, alpha, td))/2./np.power((v - tx[index])*t, 2)
    return result


def biasQestimator_shotnoise(x, paras, qshift=0, tsampling=None):
    """Qestimator model only due to the shot noise
    ; x -> Vector of binnining factors [16,32,64,128, ...]
    ; paras -> [Q0]
    ;
    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling
    return paras[0] - 1./tx


def biasQestimator_expcor(x, paras, qshift=1, tsampling=None):
    """Qestimator bias model for exponential correlation in addition to a single
        diffusing species

    x -> a size of cycle
    paras -> [Q0, td, A0, T0, alpha]

    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if not isinstance(qshift, int):
        raise TypeError("qshift is integer")
    else:
        v = qshift
    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling
    tempfcn = biasQestimator

    p1 = [paras[0], paras[1], paras[4]]
    A0  = paras[2]
    T0  = paras[3]

    result = tempfcn(x, p1, qshift=v, tsampling=t)
    result += A0*(1. + 2.*(1.- tx*t/T0 - np.exp(-tx*t/T0))/np.power(tx*t/T0, 2))

    return result


def biasQestimator_twospecies(x, paras, qshift=1, tsampling=None):
    """Qestimator bias moder for three species

    x -> Vector of binnining factors [16,32,64,128, ...]
    paras -> [Q01, td1, Q02, td2, frac, alpha]

    qshift (integer) : 1 -> Q1, 0 -> Q0
    tsampling (float) : sampling time

    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if not isinstance(qshift, int):
        raise TypeError("qshift is integer")
    else:
        v = qshift
    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling

    tempfcn = biasQestimator
    Q1, td1, Q2, td2, frac, alpha = paras[:]
    paras1=[Q1, td1, alpha]
    paras2=[Q2, td2, alpha]

    result = frac*tempfcn(tx, paras1, qshift=qshift, tsampling=tsampling)
    result += (1. - frac)*tempfcn(tx, paras2, qshift=qshift, tsampling=tsampling)

    return result


def biasQestimator_threespecies(x, paras, qshift=1, tsampling=None):
    """; x -> Vector of binnining factors [16,32,64,128, ...]
    ; paras -> [Q01, td1, Q02, td2, Q03, td3, frac1, frac2]
    ;
    ; ==> QSHIFT  (Integer) 1: calc Q1 0: calc Q0
    ; ==> alpha   (Float)  z02/w02 ;beam waist ratio
    ; ==> TSampling (Float) fundamental sampling time of FFS data in ms
    ;
    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if not isinstance(qshift, int):
        raise TypeError("qshift is integer")
    else:
        v = qshift
    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling

    tempfcn = biasQestimator
    Q1, td1, Q2, td2, Q3, td3, frac1, frac2, alpha = paras[:]
    paras1=[Q1, td1, alpha]
    paras2=[Q2, td2, alpha]
    paras3=[Q3, td3, alpha]

    result = frac1*tempfcn(tx, paras1, qshift=qshift, tsampling=tsampling)
    result += frac2*tempfcn(tx, paras2, qshift=qshift, tsampling=tsampling)
    result += (1. - frac1 - frac2)*tempfcn(tx, paras3, qshift=qshift, tsampling=tsampling)

    return result

def biasQestimator_twospecies_expcor(x, paras, qshift=1, tsampling=None):
    """Qestimator bias model for two diffusing species and a exponential
        correlation

    x -> Vector of binnining factors [16,32,64,128, ...]
    paras -> [Q01, td1, Q02, td2, frac, A0, T0, alpha]

    qshift (integer) : 1 -> Q1, 0 -> Q0
    tsampling (float) : sampling time
    """
    if isinstance(x, (int, float, list, tuple)):
        tx = np.array(x)
    elif isinstance(x, np.ndarray):
        tx = x
    else:
        raise TypeError("the data type of x is incorrect.")

    if not isinstance(qshift, int):
        raise TypeError("qshift is integer")
    else:
        v = qshift
    if tsampling is None:
        raise ValueError("tsampling is not given.")
    else:
        t = tsampling

    tempfcn = biasQestimator
    Q1, td1, Q2, td2, frac, A0, T0, alpha = paras[:]
    paras1=[Q1, td1, alpha]
    paras2=[Q2, td2, alpha]

    result = frac*tempfcn(tx, paras1, qshift=qshift, tsampling=tsampling)
    result += (1. - frac)*tempfcn(tx, paras2, qshift=qshift, tsampling=tsampling)
    result += A0*(1. + 2.*(1.- tx*t/T0 - np.exp(-tx*t/T0))/np.power(tx*t/T0, 2))

    return result
