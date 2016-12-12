import numpy as np
from mpfit.mpfit3 import mpfit

from FFSFitter import FFSFitter
import fcs_tools.fcs_models as fcs_models

class FCSFitter(FFSFitter):
    """Fit auto-correlations by fcs models for the single-color FFS time-mode
    data

    2DG : a single diffusing species in 2DG PSF
    3DG : a single diffusing species in 3DG PSF
    2DGTS : two diffusing species in 2DG PSF
    3DGTS : two diffusing species in 3DG PSF
    2DGGN : a single diffusing species in 2DG PSF in the presence of Gaussian
            Noise (exponential correlation)
    3DGGN : a single diffusing species in 3DG PSF in the presence of Gaussian
            Noise (exponential correlation)
    """
    dic_models = {"2DG":fcs_models.acf2dg, "3DG":fcs_models.acf3dg,
                    "2DGTS":fcs_models.acf2dgTS, "3DGTS":fcs_models.acf3dgTS,
                    "2DGGN":fcs_models.acf2dgGN, "3DGGN":fcs_models.acf3dgGN}
    dic_nparas = {"2DG":3, "3DG":4, "2DGTS":6, "3DGTS":7, "2DGGN":5, "3DGGN":6}

    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 4

    def __init__(self, model="2DG",
                       channels=[1],
                       paras=[0.01, 0.001, 0.],
                       fixed=[0, 0, 0]):
        if len(channels) != 1:
            raise ValueError("This FCS fitter is currently designed only for \
                single-color FFS")
        super(FCSFitter, self)._checkmodel(model)
        super(FCSFitter, self).__init__(model, channels)
        self.setModel(model, paras, fixed)

    def fit(self, t, g, gerr, tsampling=1, tseg=DATACYCLE):
        try:
            super(FCSFitter, self)._checkvalues("paras", self._paras)
            super(FCSFitter, self)._checkvalues("fixed", self._fixed)
        except:
            print("Update paras for the FCS model {}".format(self._model))
            print("Update fixed for the FCS model {}".format(self._model))
            return

        x = t
        y = g
        yerr = gerr
        if self._model in ["2DGGN", "3DGGN"]:
            fitinfo = {"tsampling":tsampling, "tseg":tseg}
        else:
            fitinfo = {"tsampling":tsampling}
            
        parinfo = [{'value':v, 'fixed':f, 'limited':[1,0], 'limits':[0.,0.]}
         					for v, f in zip(self._paras, self._fixed)]

        def myfunct(p, fjac=None, x=None, y=None, err=None, info=None):
            model =  self._fct(x, p, **info) #tsampling=info["tsampling"]
            status = 0
            return [status, (y-model)/err]

        fa = {"x":x, "y":y, "err":yerr, "info":fitinfo}
        res = mpfit(myfunct, self._paras, functkw=fa, parinfo=parinfo, maxiter=1000, \
                    quiet=1)
        yfit = self._fct(x, res.params, **fitinfo ) #tsampling=fitinfo["tsampling"]
        dof = y.size - len(self._paras) + sum(self._fixed)
        redchi2 = (np.power(y - yfit, 2)/np.power(yerr, 2)).sum()/dof
        return res, yfit, redchi2

    def setModel(self, model, paras, fixed):
        """set the FCS model, paras, and fixed together"""
        super(FCSFitter, self)._checkmodel(model)
        super(FCSFitter, self)._checkvalues("paras", paras)
        super(FCSFitter, self)._checkvalues("fixed", fixed)
        self._model = model
        self._fct = super(FCSFitter, self)._choosemodel(model)
        self._paras = paras
        self._fixed = fixed

def main():
    from readFFSfrombinfiles import readFFSfrombinfiles
    from FCSTransformer import FCSTransformer
    import matplotlib.pyplot as plt

    filenames=["A488_cal.2.001.bin"]
    data = readFFSfrombinfiles(filenames, channels=[1], frequency=100000)
    fcs = FCSTransformer()
    x = fcs.transform(data)
    fitter = FCSFitter(model="2DG", channels=[1], paras=[0.01, 0.0001, 0.], fixed=[0,0,1])

    # z0, z1, z2 = fitter.fit(x[(1,1)][0]*100000, x[(1,1)][1], x[(1,1)][2], data.sampletime)
    # z0, z1, z2 = fitter.fit(x[(0,0)], x[(1,1)][0], x[(1,1)][1], 1)
    z0, z1, z2 = fitter.fit(x[(0,0)], x[(1,1)][0], x[(1,1)][1])
    print(z0.params)
    print(z2)
    plt.plot(x[(0,0)], x[(1,1)][0])
    plt.plot(x[(0,0)], z1)
    plt.xscale("log")
    plt.show()


if __name__ == "__main__":
    main()
