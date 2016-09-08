import numpy as np
from scipy.interpolate import interp2d

"""
    t : the time
    r : (wz/w0)^2 where wz = axial beam waist and w0 = radial beam waist
    td : diffusion time 4nD/w0^2 where D = Diffusion Coefficient
"""

class Binning_Func0(object):

    def __init__(self):
        self._table = None
        self._interpolator = None

    def calculate(self, t, r, td):
        return t*0.

class Binning_Func1(object):

    def __init__(self):
        self._table = None
        self._interpolator = None

    def calculate(self, t, r, td):
        return t*1.

class Binning_Func2(object):

    def __init__(self):
        self._table = None
        self._interpolator = None

    def calculate(self, t, r, td):
        result = (4.0)/td*(r*td - np.sqrt(r*td*(r*td+t))-((td+t)            \
            *np.log(((-1.0+2.*r-2.*np.sqrt((-1.+r)*r))*(np.sqrt((-1.+r)*td) \
            + np.sqrt(r*td+t)))                                             \
            /(-np.sqrt((-1.+r)*td) + np.sqrt(r*td+t))))/(2.*np.sqrt((-1.+r)/r)))
        return td**2. * result

class Binning_Func3(object):

    def __init__(self):
        self._table = self.__loadtable()
        self._interpolator = self.__createinterpolator()

    def calculate(self, t, r, td):
        return np.power(td, 3)*t

    def __loadtable(self):
        return 0

    def __createinterpolator(self):
        interp2d(x, y, self.table)
