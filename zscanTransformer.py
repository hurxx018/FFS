import numpy as np


class zscanTransformer(object):
    """docstring for zscanTransformer.
    """
    def __init__(self, scanparas = [2.4, 0.1, 10.0515]
                     , channels = []
                     , nbins = 80
                     , start_pt = 700
                     , slice_zscans = False):
        #super(zscanTransformer, self).__init__()
        """
        scanparas = [Vpp (V), zscan_frequency (Hz)
                            , translation_distance per voltage (um/V)]
                ; trnaslation_distance is determined by the stage


        start_pt (int): the initial index where kz is kept after rebin photon
                    count data in terms of nbins
                    >> in other words, kz[start_pt:] is returned.

        slice_zscans (bool): to slice a long z-scan profile into multiple zscan
                    profiles of which each represents a single pass through the
                    sample.

        example:
            zscan = zscanTransformer()
            zscan.transform()
        """
        self._scanparas = scanparas
        self._channels = channels
        self._nbins = nbins
        self._start_pt = start_pt
        self._slice_zscans = slice_zscans
        # zscan speed : um/s
        self._zscan_speed = 2*scanparas[0]*scanparas[1]*scanparas[2]


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
        return result

    def calkz(self, data, channels, ffsfrequency):
        result = {'zscan_k':1, 'slice':self._slice_zscans}
        """
        result[0]  >> z
        result[1]  >> kz for the channel 1
        result[2]  >> kz for the channel 2 if available
        and so on.
        """
        for i in channels:
            if self._slice_zscans:
                # slice kz into mutiple pieces rather than a long kz
                temp_kz = self.sub_calkz(data[i-1])
                result[i] = self.slice_kz(temp_kz, ffsfrequency)
            else:
                result[i] = [self.sub_calkz(data[i-1])]
        # generate z-values for the position
        dist_binkz = self._zscan_speed*self._nbins/ffsfrequency
        result[0] = np.arange(result[channels[0]][0].size)*dist_binkz
        return result

    def sub_calkz(self, data):
        nk = data.size//self._nbins
        temp_data = data[:nk*self._nbins].reshape(nk, self._nbins)
        # TO DO: determine _start_pt automatically from temp_data
        # min(temp_data[:1000])
        kz = temp_data[self._start_pt:, :].sum(axis=1, dtype=np.uint16)
        return kz

    def slice_kz(self, kz, ffsfrequency):
        # time for a single pass of zscan through the sample
        time_pass = 1./self._scanparas[1]/2.
        # time for a bin in kz
        time_binkz = self._nbins*1./ffsfrequency
        # the number of bins in kz for a single pass of zscan
        nbins_pass = int(time_pass/time_binkz)
        nscans = kz.size//nbins_pass
        kz_sliced = kz[:nscans*nbins_pass].reshape(nscans, nbins_pass)
        return [kz_sliced[i] for i in range(nscans)]

    def info(self):
        return {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def __str__(self):
        temp_info = self.info()
        for key in sorted(temp_info.keys()):
            print("{0}  :  {1}".format(key[1:], temp_info[key]))
        return ""

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

    @property
    def nbins(self):
        return self._nbins

    @nbins.setter
    def nbins(self, value):
        self._nbins = value

    @property
    def start_pt(self):
        return self._start_pt

    @start_pt.setter
    def start_pt(self, value):
        self._start_pt = value

    @property
    def zscan_speed(self):
        return self._zscan_speed
