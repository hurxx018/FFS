

class FFSFitter(object):
    """a superclass for Fitters of FFS analysis

    Fits auto-correlations or cross-correlations
    """

    def __init__(self, model="", params=[], fixed=[], bounds=[], method=""):
        self.model = model
        self.params = params
        self.fixed = fixed
        self.bounds = bounds
        self.method = method


    def setModel(self, model):
        self.model = model
        return


    def setParams(self, params):
        self.params = params
        return


    def setFixed(self, fixed):
        self.fixed = fixed
        return


    def setBounds(self, bounds):
        self.bounds = bounds
        return


    def setMethod(self, method):
        self.method = method
        return


    def getModel(self):
        return self.model


    def getParams(self):
        return self.params


    def getFixed(self):
        return self.fixed


    def getBounds(self):
        return self.bounds


    def getMethod(self):
        return self.method
