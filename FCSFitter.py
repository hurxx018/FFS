import numpy as np

from . import ACF2DGDiffFitter
from . import ACF3DGDiffFitter

def FCSFitter(object):
    """Wrapper of Fitters for FCS analysis

    Fits auto-correlations or cross-correlations
    """
    def __init__(analysis="FCS", model="", params=[], fixed=[], bounds=[]):
        self._fitter = model_select(model)
        print(0)

    def fit(self, X):
        return (fitter = 0)

    def model_select(self, model):
        if model == "2DG_Diff_FCS_SINGLE":
            return ACF2DGDiffFitter.ACF2DGDiffFitter
    # @property
    # def fitmodel(self):
    #     return self._fitmodel
    #
    # @property
    # def parameters(self):
    #     return self._parameters



def main():

if __name__ == "__main__":
    main()
