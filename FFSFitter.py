

class FFSFitter(object):
    """a superclass for Fitters of FFS analysis

    Fits auto-correlations or cross-correlations
    """

    def __init__(self, model="", params=[], fixed=[], bounds=[], method=""):
        self._model = model
        self._params = params
        self._fixed = fixed
        self._bounds = bounds
        self._method = method


    def setModel(self, model):
        self._model = model
        return


    def setParams(self, params):
        self._params = params
        return


    def setFixed(self, fixed):
        self._fixed = fixed
        return


    def setBounds(self, bounds):
        self._bounds = bounds
        return


    def setMethod(self, method):
        self._method = method
        return


    def getModel(self):
        return self._model


    def getParams(self):
        return self._params


    def getFixed(self):
        return self._fixed


    def getBounds(self):
        return self._bounds


    def getMethod(self):
        return self._method
