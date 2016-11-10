import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d, interp2d

"""
    3D diffusion model and 3DG PSF model
    t : the sampling time
    r : (wz/w0)^2 where wz = axial beam waist and w0 = radial beam waist
    td : diffusion time 4nD/w0^2 where D = Diffusion Coefficient,
            n: 1 for one-photon excitation, 2 for two-photon excitation
"""

tableforB3 = "binfunctions/Binning_Func_B3_Table_For_Python"
tableforB4 = "binfunctions/Binning_Func_B4_Table_For_Python"
tableforB5 = "binfunctions/Binning_Func_B5_Table_For_Python"
tableforB6 = "binfunctions/Binning_Func_B6_Table_For_Python"
tableforB7 = "binfunctions/Binning_Func_B7_Table_For_Python"
tableforB8 = "binfunctions/Binning_Func_B8_Table_For_Python"

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
        # 3D diffusion + 3DG PSF
        result = (4.0)/td*(r*td - np.sqrt(r*td*(r*td+t))-((td+t)            \
            *np.log(((-1.0+2.*r-2.*np.sqrt((-1.+r)*r))*(np.sqrt((-1.+r)*td) \
            + np.sqrt(r*td+t)))                                             \
            /(-np.sqrt((-1.+r)*td) + np.sqrt(r*td+t))))/(2.*np.sqrt((-1.+r)/r)))
        return np.power(td, 2.) * result

class Binning_Func3(object):

    def __init__(self):
        self._table = self.loadtable()
        self._interpolator = self.createinterpolator()

    def calculate(self, t, r, td):
        dti = np.float64(0.001)
        dri = np.float64(0.01)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        ri = np.log10(r)/dri
        return np.power(td, 3)*self._interpolator(ti,ri)

    def loadtable(self):
        with open(tableforB3) as f:
            return np.fromfile(f, dtype=np.float64).reshape(8001, 401)

    def createinterpolator(self):
        return interp2d(np.arange(8001), np.arange(401)
                                    , self._table.transpose())

class Binning_Func4(object):

    def __init__(self):
        self._table = self.loadtable()
        self._interpolator = self.createinterpolator()

    def calculate(self, t, r, td):
        dti = np.float64(0.002)
        dri = np.float64(0.02)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        ri = np.log10(r)/dri
        return np.power(td, 4)*self._interpolator(ti, ri)

    def loadtable(self):
        with open(tableforB4) as f:
            return np.fromfile(f, dtype=np.float64).reshape(4001, 201)

    def createinterpolator(self):
        return interp2d(np.arange(4001), np.arange(201)
                                    , self._table.transpose())

class Binning_Func5(object):

    def __init__(self):
        self._table = self.loadtable()
        self._interpolator = self.createinterpolator()

    def calculate(self, t, r, td):
        dti = np.float64(0.01)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        if r > 100:
            rr = np.float64(100)
        else:
            rr = np.float64(r)

        if (10 < rr <= 30):
            dri = np.float64(5)
            ri = rr/dri - 2.
            return np.power(td, 5)*self._interpolator[0](ti, ri)
        elif (rr > 30 or rr > 1.1):
            dri = np.float64(10)
            ri = rr/dri
            return np.power(td, 5)*self._interpolator[1](ti, ri)
        else:
            return np.power(td, 5)*self._interpolator[2](ti)

    def loadtable(self):
        with open(tableforB5) as f:
            return np.fromfile(f, dtype=np.float64).reshape(801, 13)

    def createinterpolator(self):
        ti0 = interp2d(np.arange(801), np.arange(5)
                                , self._table[:, 1:6].transpose())
        index = np.array([0, 1, 3, 5, 6, 7, 8, 9, 10, 11, 12])
        ti1 = interp2d(np.arange(801), np.arange(11) + 0.1
                                , self._table[:, index].transpose())
        ti2 = interp1d(np.arange(801)
                                , self._table[:, 0].transpose())
        return [ti0, ti1, ti2]

# class Binning_Func6(object):
#
#     def __init__(self):
#         self._table = self.loadtable()
#         self._interpolator = self.createinterpolator()
#
#     def calculate(self, t, r, td):
#         dti = np.float64(0.002)
#         dri = np.float64(0.02)
#         ti = (np.log10(t/np.float64(td)) + 4.)/dti
#         ri = np.log10(r)/dri
#         return np.power(td, 6)*self._interpolator(ti, ri)
#
#     def loadtable(self):
#         with open(tableforB6) as f:
#             return np.fromfile(f, dtype=np.float64).reshape(801, 7)
#
#     def createinterpolator(self):
#         return interp2d(np.arange(801), np.arange(7)
#                                     , self._table.transpose())


class Binning_Func7(object):

    def __init__(self):
        self._table = self.loadtable()
        self._interpolator = self.createinterpolator()

    def calculate(self, t, r, td):
        dti = np.float64(0.01)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        return np.power(td, 7)*self._interpolator(ti)

    def loadtable(self):
        with open(tableforB7) as f:
            return np.fromfile(f, dtype=np.float64)

    def createinterpolator(self):
        return interp1d(np.arange(801), self._table.transpose())


class Binning_Func8(object):

    def __init__(self):
        self._table = self.loadtable()
        self._interpolator = self.createinterpolator()

    def calculate(self, t, r, td):
        dti = np.float64(0.01)
        ti = (np.log10(t/np.float64(td)) + 4.)/dti
        return np.power(td, 8)*self._interpolator(ti)

    def loadtable(self):
        with open(tableforB8) as f:
            return np.fromfile(f, dtype=np.float64)

    def createinterpolator(self):
        return interp1d(np.arange(801), self._table.transpose())
