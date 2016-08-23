import numpy as np
from scipy.optimize import curve_fit

from FFSFitter import FFSFitter

class FCSFitter(FFSFitter):
    """Fitter for FCS analysis

    Fits auto-correlations or cross-correlations
    """

    def __init__(self, model="2DG_SINGLE", params=None, fixed=None, \
                    bounds=None, method="lm"):
        if model == "2DG_SINGLE":
            if params != None:
                assert len(params) == 3
            if fixed != None:
                assert len(fixed) == 3
            if bounds != None:
                assert len(bounds) == 2
            FFSFitter.__init__(self, model, params, fixed, bounds, method)

        else:
            raise RuntimeError("The model {} is not available.".format(model))


    def fit(self, X, y, y_sigma=None):
        if self.getModel() == "2DG_SINGLE":
            popt, pcov = curve_fit(self.acf2DGDiff, X, y, p0=self.getParams(), \
                            sigma=y_sigma, absolute_sigma=True)#, \
                            #bounds=self._bounds)

        return popt, pcov


    def acf2DGDiff(self, x, *params):
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
        return params[0] / (1. + x / params[1]) + params[2]


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
    print(acf.getInfo())
    xx = acf.transform(ffsdata)

    acffit = FCSFitter(model="2DG_SINGLE", params=[0.001, 0.0001, 0.00001])
    a0, a1 = acffit.fit(xx[(1,1)].time, xx[(1,1)].correlations, xx[(1,1)].correlations_stderr)

    print(a0)
    print("-------")
    print(a1)

if __name__ == "__main__":
    main()
