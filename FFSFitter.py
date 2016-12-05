
class FFSFitter(object):
    """a superclass for Fitters of FFS analysis

    """
    dic_models = {}
    dic_nparas = {}

    def __init__(self, model="", channels=[], paras=[], fixed=[]):
        self._model = model
        self._channels = channels
        self._paras = paras
        self._fixed = fixed
        self._fct = None

    def getmodelfct(self):
        """get a model function"""
        return self._fct

    def _choosemodel(self, model):
        return self.dic_models[model]

    def _checkmodel(self, model):
        if model not in self.dic_models:
            for key in sorted(self.dic_models):
                print("Use {} for {}".format(key, self.dic_models[key].__name__))
            print("")
            raise ValueError("{} model is not available.".format(model))

    def _checkvalues(self, name, value):
        if not isinstance(value, (list, tuple)):
            raise TypeError("the type of {} is neither a list or a tuple".format(name))
        np = self.dic_nparas[self._model]
        if len(value) != np:
            raise ValueError("The number of elements in {} should be {} \
            for the MSQ model {}.".format(name, np, self._model))

    def info(self):
        return {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def __str__(self):
        for key, value in self.info().items():
            print("{0}  :  {1}".format(key[1:], value))
        return ""

    @property
    def model(self):
        return self._model
    @model.setter
    def model(self, value):
        self._checkmodel(value)
        self._model = value
        self._fct = self._choosemodel(model)

    @property
    def channels(self):
        """Fit method for Fitter"""
        return self._channels

    @channels.setter
    def channels(self, value):
        if isinstance(value, (list, tuple)):
            self._channels = value
        else:
            raise TypeError("the type of channels is neither a list or a tuple")

    @property
    def paras(self):
        return self._paras
    @paras.setter
    def paras(self, value):
        self._checkvalues("paras", value)
        self._paras = value

    @property
    def fixed(self):
        return self._fixed
    @fixed.setter
    def fixed(self, value):
        self._checkvalues("fixed", value)
        self._fixed = value
