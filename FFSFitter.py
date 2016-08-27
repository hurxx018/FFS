
class FFSFitter(object):
    """a superclass for Fitters of FFS analysis

    inputs:
    analysis = "FCS", "TIFCA", "TSFCA" and so on

    model:
    (1) FCS analysis >> "2DG_FCS_SINGLE"
    (2) TIFCA analysis >>
    (3) TSFCA analysis >>

    """

    def __init__(self, analysis="", model="", method=""):
        self._analysis = analysis
        self._model = model
        self._method = method


    @property
    def analysis(self):
        """A type of FFS anaysis"""
        return self._analysis

    @analysis.setter
    def analysis(self, value):
        self._analysis = value


    @property
    def model(self):
        """Fit model for Fitter"""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value


    @property
    def method(self):
        """Fit method for Fitter"""
        return self._method

    @method.setter
    def method(self, value):
        self._method = value


    # def setModel(self, model):
    #     self.model = model
    #     return
    #
    #
    # def setParams(self, params):
    #     self.params = params
    #     return
    #
    #
    # def setFixed(self, fixed):
    #     self.fixed = fixed
    #     return
    #
    #
    # def setBounds(self, bounds):
    #     self.bounds = bounds
    #     return
    #
    #
    # def setMethod(self, method):
    #     self.method = method
    #     return
    #
    #
    # def getModel(self):
    #     return self.model
    #
    #
    # def getParams(self):
    #     return self.params
    #
    #
    # def getFixed(self):
    #     return self.fixed
    #
    #
    # def getBounds(self):
    #     return self.bounds
    #
    #
    # def getMethod(self):
    #     return self.method
