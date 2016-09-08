

class FFSTransformer(object):
    """Superclass of transformers in FFS analysis


    """

    def __init__(self, segmentlength=32768*4
                     , channels = []):
        self._segmentlength = segmentlength
        self._channels = channels

    def fit(self):
        return self

    def transform(self):
        """this is defined in each subclass"""
        pass

    def info(self):
        return {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def __str__(self):
        for key, value in self.info().items():
            print("{0}  :  {1}".format(key[1:], value))
        return ""

    @property
    def segmentlength(self):
        return self._segmentlength

    @segmentlength.setter
    def segmentlength(self, value):
        self._segmentlength = value

    @property
    def channels(self):
        return self._channels

    @channels.setter
    def channels(self, channels):
        self._channels = channels
