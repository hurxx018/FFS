import numpy as np

from zscan_tools.zscanPSFGLsp import volPSFGLsp, volPSFsqrGLsp

def g2ratio_SLAB(x, paras):
    """Calculate Qs in the SLAB model with the mGL PSF
    x >> thickness of slab
    paras = [zR, y]

    """
    p1 = paras
    g2inf = volPSFsqrGLsp(p1, None, None, geo=-1)/volPSFGLsp(p1, None, None, geo=-1)
    g2fin = volPSFsqrGLsp(p1, -x/2., x/2., geo=0) \
                /volPSFGLsp(p1, -x/2., x/2., geo=0)
    return g2fin/g2inf


class zscan_g2ratio(object):
    """docstring for zscan_g2ratio.
    (L, Q) arrays

    filename >> (L, Q) text file
    paras >> [eta, z0]

    Get the result
    result()
    save() >> save the result into a csv file

    """
    def __init__(self, filename, paras):
        super(zscan_g2ratio, self).__init__()
        self._filename = filename
        self._data = np.loadtxt(self._filename, dtype=np.float32)

        x0, x1 = self.g2corrected(paras)

        self._g2ratio = x0
        self._qcorrected = x1

    def g2corrected(self, paras):
        x0 = g2ratio_SLAB(self._data[:, 0], paras)
        x1 = self._data[:, 1]/x0
        return x0, x1

    def result(self):
        return np.hstack((self._data, self._g2ratio[:, np.newaxis],
                                    self._qcorrected[:, np.newaxis]))

    def save(self, filename=None):
        res = self.result()
        if filename:
            np.savetxt(filename, res, fmt="%5.4f, %5.4f, %5.4f, %5.4f", newline='\n')
        else:
            np.savetxt("LQ_g2corrected.csv", res, fmt="%5.4f, %5.4f, %5.4f, %5.4f", newline='\n')

    @property
    def g2ratio(self):
        return self._g2ratio

    @property
    def qcorrected(self):
        return self._qcorrected
