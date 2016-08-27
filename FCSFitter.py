import numpy as np
#from scipy.optimize import curve_fit
from lmfit import Parameters, Model

from FFSFitter import FFSFitter

class FCSFitter(FFSFitter):
    """Fitter for FCS analysis

    Fits auto-correlations or cross-correlations
    """

    def __init__(self, analysis="FCS"
                     , model="2DG_FCS_SINGLE"
                     , params=[0.001, 0.0001, 0.]
                     , fixed=[0, 0, 0]
                     , bounds=[[0, np.inf], [0, np.inf], [0, np.inf]]
                     , method="leastsq"):

        if analysis=="FCS" and model == "2DG_FCS_SINGLE":
            assert len(params) == 3
            assert len(fixed) == 3
            assert len(bounds) == 3
            for bound in bounds:
                assert len(bound) == 2
            FFSFitter.__init__(self, analysis, model, method)
            self._fitmodel = Model(self.acf2DGDiff
                                    , independent_vars=['x']
                                    , param_names=["g0", "td", "offset"])
            #
            # self._parameters = Parameters()
            for i, key in enumerate(self._fitmodel.param_names): #["g0", "td", "offset"]
                self._fitmodel.set_param_hint(key, value=params[i]
                                                , vary= (not bool(fixed[i]))
                                                , min=bounds[i][0]
                                                , max=bounds[i][1])
            self._parameters = self._fitmodel.make_params()


        else:
            raise RuntimeError("The analysis and the model should be {} and {}."
                                .format(["FCS", "2DG_FCS_SINGLE"]))

    # TO DO: Employ the LMFIT package
    def fit(self, X, y, y_sigma=None):
        if y_sigma is not None:
            return self.fitmodel.fit(y, params=self._parameters
                                      , weights=1./y_sigma
                                      , method=self.method
                                      , x=X)
        else:
            return self.fitmodel.fit(y, params=self._parameters
                                      , method=self.method
                                      , x=X)


        # if self.getModel() == "2DG_FCS_SINGLE":
        #     popt, pcov = curve_fit(self.acf2DGDiff, X, y, p0=self.getParams(), \
        #                     sigma=y_sigma, absolute_sigma=True)#, \
        #                     #bounds=self._bounds)
        #
        # return popt, pcov


#    def acf2DGDiff(self, x, *params):
    def acf2DGDiff(self, x, g0, td, offset):
        """Auto-correlation function for a single species
        modeled with simple diffusion and 2D Gaussian Point-Spred Function (PSF)

        Parameters
        ----------
        x : array
        params : [g0, td, offset]

        Returns
        -------
        array
        """
#        return params[0] / (1. + x / params[1]) + params[2]
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

    acffit = FCSFitter(model="2DG_SINGLE", params=[0.001, 0.0001, 0.00001])
    a0, a1 = acffit.fit(xx[(1,1)].time, xx[(1,1)].correlations, xx[(1,1)].correlations_stderr)

    print(a0)
    print("-------")
    print(a1)

if __name__ == "__main__":
    main()
