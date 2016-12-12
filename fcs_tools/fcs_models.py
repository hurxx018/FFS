from numpy import ndarray, sqrt, exp, power, zeros, ones


"""Module for models for Single color FCS analysis
Function Names:
    A single diffusing species in 2DG PSF  >> acf2dg
    A single diffusing species in 3DG PSF  >> acf3dg
    Two diffusing species in 2DG PSF  >> acf2dgTS
    Two diffusing species in 3DG PSF  >> acf3dgTS
    A single diffusing species in 2DG PSF in the presence of Gaussian Noise
    (exponential correlation) >> acf2dgGN
    A single diffusing species in 3DG PSF in the presence of Gaussian Noise
    (exponential correlation) >> acf3dgGN
"""

CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
DATACYCLE = CYCLE * 4


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

def acf2dgTS(x, paras, tsampling=1):
    """Autocorrelation function of two diffusing species in 2DG PSF

    paras : [G0_1, td_1, G0_2, td_2, f1, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 6:
        raise ValueError("The number of elements in paras should be 3")
    return paras[4]*paras[0]/(1. + t/paras[1]) \
            + (1.-paras[4])*paras[2]/(1. + t/paras[3]) \
            + paras[5]

def acf3dgTS(x, paras, tsampling=1):
    """Autocorrelation function of a single diffusing species in 3DG PSF

    paras : [G0_1, td_1, G0_2, td_2, f1, r^2, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 7:
        raise ValueError("The number of elements in paras should be 4")
    return paras[4]*paras[0]/(1. + t/paras[1])/sqrt(1. + t/paras[1]/paras[5]) \
        + (1. - paras[4])*paras[2]/(1. + t/paras[3])/sqrt(1. + t/paras[3]/paras[5]) \
        + paras[6]

def b2gn(x, paras, tsampling=1):
    """2nd order binning function for the Gaussian Noise(exponential correlation)"""
    if tsampling == 1:
        t = _checkx(x)
    else:
        t = _checkx(x)*tsampling
    if len(paras) != 1:
        raise ValueError("The number of elements in paras should be 1")
    return 2*power(paras[0], 2)*(-1. + t/paras[0] + exp(-t/paras[0]))

def acf2dgGN(x, paras, tsampling=1, tseg=DATACYCLE):
    """Autocorrelation function of a single diffusing species in 2DG PSF in the
    presence of exponential correlation (Gaussian Noise) and its estimation bias
    correction

    paras : [G0, td, A0, T0, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
        tseg_ = ones(t.size) * tseg
    else:
        t = _checkx(x)*tsampling
        tseg_ = ones(t.size) * tseg * tsampling

    if len(paras) != 5:
        raise ValueError("The number of elements in paras should be 5")
    result = paras[0]/(1. + t/paras[1]) + paras[2]*exp(-t/paras[3]) + paras[4]

    fcn = b2gn

    indices = (t <= tseg_/2)
    temp_tseg_ = tseg_[indices]
    temp_t_ = t[indices]
    result[indices] += paras[2]*(fcn(temp_tseg_, [paras[3]], tsampling)
            + fcn((temp_tseg_ - 2*temp_t_ ), [paras[3]], tsampling)
            - 2.*fcn(temp_t_ , [paras[3]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)

    indices = (t > tseg_/2)
    if indices.any():
        temp_tseg_ = tseg_[indices]
        temp_t_ = t[indices]
        result[indices] += paras[2]*(fcn(temp_tseg_, [paras[3]], tsampling)
            + fcn((2*temp_t_ - temp_tseg_), [paras[3]], tsampling)
            - 2.*fcn(temp_t_ , [paras[3]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)

    return result

def acf3dgGN(x, paras, tsampling=1, tseg=DATACYCLE):
    """Autocorrelation function of a single diffusing species in 3DG PSF in the
    presence of exponential correlation (Gaussian Noise) and its estimation bias
    correction

    paras : [G0, td, A0, T0, r^2, offset]
    """
    if tsampling == 1:
        t = _checkx(x)
        tseg_ = ones(t.size) * tseg
    else:
        t = _checkx(x)*tsampling
        tseg_ = ones(t.size) * tseg * tsampling

    if len(paras) != 6:
        raise ValueError("The number of elements in paras should be 6")
    result = paras[0]/(1. + t/paras[1])/sqrt(1. + t/paras[1]/paras[4]) \
            + paras[2]*exp(-t/paras[3]) + paras[5]
# paras[0]/(1. + t/paras[1])/sqrt(1. + t/paras[1]/paras[2]) + paras[3]
    fcn = b2gn

    indices = (t <= tseg_/2)
    temp_tseg_ = tseg_[indices]
    temp_t_ = t[indices]
    result[indices] += paras[2]*(fcn(temp_tseg_, [paras[3]], tsampling)
            + fcn((temp_tseg_ - 2*temp_t_ ), [paras[3]], tsampling)
            - 2.*fcn(temp_t_ , [paras[3]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)

    indices = (t > tseg_/2)
    if indices.any():
        temp_tseg_ = tseg_[indices]
        temp_t_ = t[indices]
        result[indices] += paras[2]*(fcn(temp_tseg_, [paras[3]], tsampling)
            + fcn((2*temp_t_ - temp_tseg_), [paras[3]], tsampling)
            - 2.*fcn(temp_t_ , [paras[3]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)

    return result

def estbias_gn(x, paras, tsampling=1, tseg=DATACYCLE):
    """Estimation Bias of autocorrelation due to Gaussian Noise

    paras : [A0, T0]
    """
    if tsampling == 1:
        t = _checkx(x)
        tseg_ = ones(t.size) * tseg
    else:
        t = _checkx(x)*tsampling
        tseg_ = ones(t.size) * tseg * tsampling

    if len(paras) != 2:
        raise ValueError("The number of elements in paras should be 2")

    fcn = b2gn
    result = zeros(t.size)

    indices = (t <= tseg_/2)
    temp_tseg_ = tseg_[indices]
    temp_t_ = t[indices]
    result[indices] += paras[0]*(fcn(temp_tseg_, [paras[1]], tsampling)
            + fcn((temp_tseg_ - 2*temp_t_ ), [paras[1]], tsampling)
            - 2.*fcn(temp_t_ , [paras[1]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)

    indices = (t > tseg_/2)
    if indices.any():
        temp_tseg_ = tseg_[indices]
        temp_t_ = t[indices]
        result[indices] += paras[0]*(fcn(temp_tseg_, [paras[1]], tsampling)
                + fcn((2*temp_t_ - temp_tseg_), [paras[1]], tsampling)
                - 2.*fcn(temp_t_ , [paras[1]], tsampling))/2./power(temp_tseg_ - temp_t_, 2)
    return result



def main():
    import matplotlib.pyplot as plt
    from numpy import arange
    x = (arange(110000) + 1)
    plt.plot(x*0.00001, acf2dg(x, [1, 0.001, 0.], 0.00001), label='2DG')
    plt.plot(x*0.00001, acf3dg(x, [1, 0.001, 25., 0.], 0.00001), label='3DG')
    plt.plot(x*0.00001, acf2dgGN(x, [1, 0.001, 0.1, 1., 0.], 0.00001, 32768*4), label='2DGGN')
    plt.plot(x*0.00001, acf3dgGN(x, [1, 0.001, 0.1, 1., 25., 0.], 0.00001, 32768*4), label='3DGGN')
    plt.plot(x*0.00001, acf2dgTS(x, [1, 0.01, 1, 0.001, 0.5, 0.], 0.00001), label='2DGTS')
    plt.plot(x*0.00001, acf3dgTS(x, [1, 0.01, 1, 0.001, 0.5, 25., 0.], 0.00001), label='3DGTS')
    plt.xscale("log")
    plt.legend()
    plt.show()

if __name__=="__main__":
    main()
