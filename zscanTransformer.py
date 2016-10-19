import numpy as np


class zscanTransformer(object):
    """docstring for zscanTransformer.
    """
    def __init__(self, scanparas = [2.4, 0.1, 10.0515]
                     , channels = []
                     , nbins = 80):
        #super(zscanTransformer, self).__init__()
        """
        scanparas = [Vpp (V), zscan_frequency (Hz)
                            , translation_distance per voltage (um/V)]
                ; trnaslation_distance is determined by the stage
        """
        self._scanparas = scanparas
        self._channels = channels
        self._nbins = nbins


    def transform(self, X):
        if self._channels == []:
            data = X.data
            channels = X.channels
            ffsfrequency = X.frequency
            result = self.calkz(data, channels, ffsfrequency)
        else:
            data = X.data
            ffsfrequency = X.frequency
            result = self.calkz(data, self._channels, ffsfrequency)

    def calkz(self, data, channels, ffsfrequency):
        result = {}
        for i in channels:
            result[i] = self.sub_calkz(data[i-1], ffsfrequency)
        return result

    def sub_calkz(self, data, ffsfrequency):
        nk = data.size//self._nbins
        temp_data = data[:nk*self._nbins]
        kz = temp_data.reshape(nk, self._nbins).sum(axis=1, dtype=uint16)

        return kz



    @property
    def scanparas(self):
        return self._scanparas


    @property
    def channels(self):
        return self._channels

    @channels.setter
    def channels(self, channels):
        if isinstance(channels, list):
            self._channels = channels
        else:
            raise ValueError("The data type of channels is list.")
