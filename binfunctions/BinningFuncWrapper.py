import sys, os
import numpy as np

from scipy.interpolate import interp2d

tableforB3 = "binfunctions/Binning_Func_B3_Table_For_Python"
tableforB4 = "binfunctions/Binning_Func_B4_Table_For_Python"
tableforB5 = "binfunctions/Binning_Func_B5_Table_For_Python"
tableforB6 = "binfunctions/Binning_Func_B6_Table_For_Python"
tableforB7 = "binfunctions/Binning_Func_B7_Table_For_Python"
tableforB8 = "binfunctions/Binning_Func_B8_Table_For_Python"


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
        if order not in range(4):
            raise ValueError("{} order is not available".format(order))
        self._order = order
        self._table = self.loadtable(order)
        self._interpolator = self.createinterpolator(order)
        self._fcn = self.function(order)

    def calculate(self, t, r, td):
        return self._fcn(t, r, td)

    def function(self, order):
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
        return np.power(td, 2) * result

    def Binning_Func_B3(self, t, r, td):
        dti = np.float64(0.001)
        dri = np.float64(0.01)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        ri = np.log10(r)/dri
        return np.power(td, 3)*self._interpolator(ti,ri)

    def Binning_Func_B4(self, t, r, td):
        dti = np.float64(0.002)
        dri = np.float64(0.02)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        ri = np.log10(r)/dri
        return np.power(td, 4)*self._interpolator(ti, ri)


    def loadtable(self, order):

        if order in [0, 1, 2]:
            return None
        elif order == 3:
            #TO DO: how to read Binning_Func_B3_Table_For_Python
            #       which is in the same folder of BinningFuncWrapper
            with open(tableforB3, "r") as f:
                return np.fromfile(f, dtype=np.float64).reshape(8001, 401)
        elif order == 4:
            with open(tableforB4) as f:
                return np.fromfile(f, dtype=np.float64).reshape(4001, 201)

    def createinterpolator(self, order):
        if order in [0, 1, 2]:
            return None
        elif order == 3:
            return interp2d(np.arange(8001), np.arange(401)
                                    , self._table.transpose())
        elif order == 4:
            return interp2d(np.arange(4001), np.arange(201)
                                    , self._table.transpose())

    @property
    def fcn(self):
        return self._fcn

    @fcn.setter
    def fcn(self, order):
        self._fcn = self.function(order)

    @property
    def table(self):
        return self._table
