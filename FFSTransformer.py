

class FFSTransformer(object):
    """Superclass of transformers in FFS analysis

    Inputs:
        segmentlength >> lenghth per segment
        channels      >> a list of channel labels e.g. [1, 2]



    """

    def __init__(self, segmentlength=32768*4
                     , channels = []):
        self.segmentlength = segmentlength
        self.channels = channels

    # def fit(self):
    #     return self

    def transform(self):
        """ An instacne method of FFS estimators from tmdata.
        >> this method is defined in each subclass.

        """
        pass

    def info(self):
        return {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def __str__(self):
        for key, value in self.info().items():
            print("{0:30s}  :  {1}".format(key, value))

    @property
    def segmentlength(self):
        return self.__segmentlength
    @segmentlength.setter
    def segmentlength(self, value):
        assert isinstance(value, (float, int)), \
            "TypeError: segmentlength should be either integer or float."
        self.__segmentlength = value

    @property
    def channels(self):
        return self.__channels
    @channels.setter
    def channels(self, channels):
        assert isinstance(channels, (list, tuple)), \
            "Type of channels should be either list or tuple."
        self.__channels = list(channels)


    # Auxillary methods
    @property
    def seglen(self):
        return self.segmentlength
