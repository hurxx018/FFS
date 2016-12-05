import numpy as np
from mpfit.mpfit3 import mpfit

from FFSFitter import FFSFitter
import fcs_tools.fcs_models as fcs_models

class FCSFitter(FFSFitter):
    """Fit auto-correlations by fcs models for the single-color FFS time-mode
    data

    """
    dic_models = {"2DG":0, "3DG":2}
    dic_nparas = {"2DG":3, "3DG":4}

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

    def fit(self, t, g, gerr, tsampling):
        try:
            super(FCSFitter, self)._checkvalues("paras", self._paras)
            super(FCSFitter, self)._checkvalues("fixed", self._fixed)
        except:
            print("Update paras for the MSQ model {}".format(self._model))
            print("Update fixed for the MSQ model {}".format(self._model))
            return

        x = t
        y = g
        yerr = gerr

        fitinfo = {"tsampling":tsampling }
        parinfo = [{'value':v, 'fixed':f, 'limited':[1,0], 'limits':[0.,0.]}
         					for v, f in zip(self._paras, self._fixed)]

        def myfunct(p, fjac=None, x=None, y=None, err=None, info=None):
            model =  self._fct(x, p, tsampling=info["tsampling"])
            status = 0
            return [status, (y-model)/err]

        fa = {"x":x, "y":y, "err":yerr, "info":fitinfo}
        # res = mpfit(myfunct, self._paras, functkw=fa, parinfo=parinfo, maxiter=1000, \
        #             quiet=1)
        yfit = self._fct(x, res.params, tsampling=fitinfo["tsampling"])
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
    a = FCSFitter()
    print(a)
    b = FCSFitter('3DG', [1], [0.01, 0.001, 25., 0], [0,0,0,0])
    print(b)

if __name__ == "__main__":
    main()
