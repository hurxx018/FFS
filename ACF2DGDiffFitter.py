import numpy as np
from lmfit import Parameters, Model

from FFSFitter import FFSFitter

class ACF2DGDiffFitter(FFSFitter):
    """Fitter for FCS analysis of 2DG PSF, Diffusion, and single species

    Fits auto-correlations
    """

    def __init__(self, analysis="FCS"
                     , model="2DG_Diff_FCS_SINGLE"
                     , params=[0.001, 0.0001, 0.]
                     , fixed=[0, 0, 0]
                     , bounds=[[0, np.inf], [0, np.inf], [0, np.inf]]
                     , method="leastsq"):

        if analysis == "FCS" and model == "2DG_Diff_FCS_SINGLE":
            assert len(params) == 3
            assert len(fixed) == 3
            assert len(bounds) == 3
            for bound in bounds:
                assert len(bound) == 2
            FFSFitter.__init__(self, analysis, model, method)
            self._fitmodel = Model(self.acf2DGDiff
                                    , independent_vars=['x']
                                    , param_names=["g0", "td", "offset"])
            for i, key in enumerate(self._fitmodel.param_names):
                self._fitmodel.set_param_hint(key, value=params[i]
                                                , vary= (not bool(fixed[i]))
                                                , min=bounds[i][0]
                                                , max=bounds[i][1])
            self._parameters = self._fitmodel.make_params()


        else:
            raise RuntimeError("The analysis and the model should be {} and {}."
                                .format(["FCS", "2DG_Diff_FCS_SINGLE"]))


    def fit(self, X, y, y_sigma=None, fit_kws=None):
        """Fit the model to the data.

        Parameters
        ----------
        X : array-like time for FCS analysis
        y : array-like Auto-correlation
        y_sigma : array-like standard-deviation

        fit_kws: dict
            default fitting options, such as xtol and maxfev, for scipy.optimizer
            fit_kws = {"xtol" : 0.0000001}


        Returns
        -------
        lmfit.ModelResult


        """
        if y_sigma is not None:
            return self.fitmodel.fit(y, params=self.parameters
                                      , weights=1./y_sigma
                                      , method=self.method
                                      , fit_kws=fit_kws
                                      , x=X)
        else:
            return self.fitmodel.fit(y, params=self.parameters
                                      , method=self.method
                                      , fit_kws=fit_kws
                                      , x=X)


    def acf2DGDiff(self, x, g0, td, offset):
        """Auto-correlation function for a single species
        modeled with simple diffusion and 2D Gaussian Point-Spred Function (PSF)

        Parameters
        ----------
        x : array
        params : g0, td, offset

        Returns
        -------
        array
        """
        return g0 / (1. + x / td) + offset


    @property
    def fitmodel(self):
        return self._fitmodel

    @property
    def parameters(self):
        return self._parameters



def main():
    import readFFSfrombinfiles as rffs
    import calfcsTransformer as fcs

    filename2 = "A488_cal.2.001.bin"

    ffsdata = rffs.readFFSfrombinfiles([filename2], [1], frequency=100000)
    print("filenames : ")
    for x in ffsdata.getFilenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.getChannels())

    acf = fcs.calfcsTransformer(channels=[1])
    print(acf)
    xx = acf.transform(ffsdata)

    acffit = ACF2DGDiffFitter(analysis="FCS"
                        , model="2DG_Diff_FCS_SINGLE"
                        , params=[0.001, 0.0001, 0.00001])
    result = acffit.fit(xx[(1,1)].time, xx[(1,1)].correlations
                            , y_sigma= xx[(1,1)].correlations_stderr
                            , fit_kws={'xtol':0.00001})

    print(result.fit_report())
    # print("-------")
    # print(a1)

if __name__ == "__main__":
    main()
