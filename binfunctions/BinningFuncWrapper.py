import sys, os
import numpy as np

from scipy.interpolate import interp2d




class BinningFuncWrapper(object):
    """Wrapper for the Binnning Functions

    order : the order of Binning Functions
    t : the time
    r : (wz/w0)^2 where wz = axial beam waist and w0 = radial beam waist
    td : diffusion time 4nD/w0^2 where D = Diffusion Coefficient

    Example
    -------
    bf2 = BinningFuncWrapper(2)
    bf2.calculate(1., 25., 1)
    """

    # TO DO: use
    #self.directory = "C:\\Users\\MLab\\Google Drive\\Python Application\FFS\\binfunctions\\"

    def __init__(self, order):
        self._order = order
        self._table = self.__loadtable(order)
        self._interpolator = self.__createinterpolator(order)
        self._fcn = self.__function(order)

    def calculate(self, t, r, td):
        return self.fcn(t, r, td)

    def __function(self, order):
        if order == 0:
            return self.Binning_Func_B0
        elif order == 1:
            return self.Binning_Func_B1
        elif order == 2:
            return self.Binning_Func_B2
        elif order == 3:
            return self.Binning_Func_B3


    def Binning_Func_B0(self, t, r, td):
        return t*0.

    def Binning_Func_B1(self, t, r, td):
        return t*1.

    def Binning_Func_B2(self, t, r, td):
        result = (4.0)/td*(r*td - np.sqrt(r*td*(r*td+t))-((td+t)            \
                *np.log(((-1.0+2.*r-2.*np.sqrt((-1.+r)*r))*(np.sqrt((-1.+r)*td) \
                +np.sqrt(r*td+t)))                                          \
                /(-np.sqrt((-1.+r)*td) + np.sqrt(r*td+t))))/(2.*np.sqrt((-1.+r)/r)))
        return td**2. * result
    # TO DO
    def Binning_Func_B3(self, t, r, td):
        # = 0.001
        # = 0.01
        #
        result = 0
        return td**3. * result


    def __loadtable(self, order):

        if order in [0, 1, 2]:
            return None
        elif order == 3:
            #TO DO: how to read Binning_Func_B3_Table_For_Python
            #       which is in the same folder of BinningFuncWrapper
            temp_table = np.fromfile("C:\\Users\\MLab\\Google Drive\\Python Application\FFS\\binfunctions\\Binning_Func_B3_Table_For_Python")

            return temp_table



    @property
    def fcn(self):
        return self._fcn

    @fcn.setter
    def fcn(self, order):
        self._fcn = self.__function(order)

    @property
    def table(self):
        return self._table
