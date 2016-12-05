import numpy as np
from mpfit.mpfit3 import mpfit

import msq_tools.msq_models as msq_models
from FFSFitter import FFSFitter

class MSQFitter(FFSFitter):
    """Fit MSQ estimators by MSQ models for single-color FFS time-mode data

    MSQ models:
    SS : 3DG PSF + a single diffusing species
        => paras = [Q0, td, alpha]

    TS : 3DG PSF + two diffusing species
        => paras = [Q0_1, td_1, Q0_2, td_2, frac1, alpha]

    THS : 3DG PSF + three diffusing species
        => paras = [Q0_1, td_1, Q0_2, td_2, Q0_3, td_3, frac1, frac2, alpha]

    SN : only due to the shot noise
        => paras = [Q0]

    SSEC : 3DG + a single diffusing species + exponential correlation
        => paras = [Q0, td, A0, T0, alpha]

    TSEC : 3DG + two diffusion species + exponential correlation
        => paras = [Q0_1, td_1, Q0_2, td_2, frac1, A0, T0, alpha]

    """

    dic_models = {"SS":msq_models.biasQestimator,
                    "TS":msq_models.biasQestimator_twospecies,
                    "THS":msq_models.biasQestimator_threespecies,
                    "SN":msq_models.biasQestimator_shotnoise,
                    "SSEC":msq_models.biasQestimator_expcor,
                    "TSEC":msq_models.biasQestimator_twospecies_expcor}

    dic_nparas = {"SS":3, "TS":6, "THS":9, "SN":1, "SSEC":5, "TSEC":8}

    def __init__(self, model="SS", \
                       channels=[1], \
                       paras=[0.01, 0.001, 25.], \
                       fixed=[0, 0, 0]): #, frequency=None, qshift=1):
        if len(channels) != 1:
            raise ValueError("This MSQ fitter is currently designed only for \
                single-color FFS")
        super(MSQFitter, self)._checkmodel(model)
        super(MSQFitter, self).__init__(model, channels)
        self.setModel(model, paras, fixed)

    def fit(self, t, qest, qesterr, qshift, tsampling):
        try:
            super(MSQFitter, self)._checkvalues("paras", self._paras)
            super(MSQFitter, self)._checkvalues("fixed", self._fixed)
        except:
            print("Update paras for the MSQ model {}".format(self._model))
            print("Update fixed for the MSQ model {}".format(self._model))
            return

        x = t
        y = qest
        yerr = qesterr

        fitinfo = {"qshift":qshift, "tsampling":tsampling }
        parinfo = [{'value':v, 'fixed':f, 'limited':[1,0], 'limits':[0.,0.]}
         					for v, f in zip(self._paras, self._fixed)]

        def myfunct(p, fjac=None, x=None, y=None, err=None, info=None):
            model =  self._fct(x, p, qshift=info["qshift"], tsampling=info["tsampling"])
            status = 0
            return [status, (y-model)/err]

        fa = {"x":x, "y":y, "err":yerr, "info":fitinfo}
        res = mpfit(myfunct, self._paras, functkw=fa, parinfo=parinfo, maxiter=1000, \
                    quiet=1)
        yfit = self._fct(x, res.params, qshift=fitinfo["qshift"], tsampling=fitinfo["tsampling"])
        dof = y.size - len(self._paras) + sum(self._fixed)
        redchi2 = (np.power(y - yfit, 2)/np.power(yerr, 2)).sum()/dof
        return res, yfit, redchi2

    def setModel(self, model, paras, fixed):
        """set the MSQ model, paras, and fixed together"""
        super(MSQFitter, self)._checkmodel(model)
        super(MSQFitter, self)._checkvalues("paras", paras)
        super(MSQFitter, self)._checkvalues("fixed", fixed)
        self._model = model
        self._fct = super(MSQFitter, self)._choosemodel(model)
        self._paras = paras
        self._fixed = fixed


def main():
    from readFFSfrombinfiles import readFFSfrombinfiles
    from MSQTransformer import MSQTransformer
    import matplotlib.pyplot as plt

    filenames=["A488_cal.2.001.bin"]
    data = readFFSfrombinfiles(filenames, channels=[1], frequency=100000)
    msq = MSQTransformer(segmentlength=32768*8*4, binfactors=(2**np.arange(13)).tolist(), tshift=1)
    x = msq.transform(data)
    fitter = MSQFitter(model="SS", channels=[1], paras=[0.01, 0.0001, 25.], fixed=[0,0,1])


    print(x[1][0])
    print(x[1][1])
    print(x[1][2])
    z0, z1, z2 = fitter.fit(x[1][0]*100000., x[1][1], x[1][2], 1,  data.sampletime)
    print(z2)
    print(z0.params)
    plt.plot(x[1][0], x[1][1])
    plt.plot(x[1][0], z1)
    plt.xscale("log")
    plt.show()
if __name__=="__main__":
    main()
