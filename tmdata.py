

class tmdata(object):
    """
    A superclass for the FFS data
    input variables:
    filenames : a python list of filenames
    channels : a python list of channels in FFS data.
    frequency: an integer of frequency of FFS data collection

    Instance Methods:
    Setters: setFilenames, setChannels, setFrequency
    Getters: getFilenames, getChannels, getNChannels, getFrequency

    """

    def __init__(self, filenames=[], channels=[], frequency=1):
        if isinstance(filenames, list):
            self._filenames = filenames
        else:
            raise TypeError("filenames is not a python list.")

        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")

        if isinstance(frequency, int):
            self._frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")
        self._data = [None]*len(channels)

    def rebin(self, cha=1, lbin=1):
        self._checkchalbin(cha, lbin)
        if lbin == 1:
            return self._data[cha-1], self._frequency
        else:
            bindatasize = self._data[cha-1].size//lbin
            k = (self._data[cha-1][:bindatasize*lbin]
                     .reshape(bindatasize, lbin).sum(axis=1))
            return k, self._frequency/lbin

    def _checkdata(self):
        if self.nchannels > 1:
            for i in range(self.nchannels - 1):
                assert len(self._data[i]) == len(self._data[i+1]), \
                    "The number of data points does not match."

    def _checkchalbin(self, cha, lbin):
        if not (isinstance(cha, int) and (cha in range(1, self.nchannels+1))):
            raise ValueError("{} is not available".format(cha))
        if not isinstance(lbin, int):
            raise TypeError('The type of lbin is int.')

    @property
    def filenames(self):
        pass # defined in subclasses

    @filenames.setter
    def filenames(self, filenames):
        pass # defined in subclasses

    @property
    def channels(self):
        return self._channels

    @channels.setter
    def channels(self, channels):
        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")

    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, frequency):
        if isinstance(frequency, int):
            self._frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, data):
        if not isinstance(data, list):
            raise TypeError("The type of data is list.")
        if len(data) != self.nchannels:
            print("The number of channels is required to be {}".format(self.nchannels))
            raise ValueError("The number of channels is not correct")
        self._data = data

    @property
    def nchannels(self):
        return len(self._channels)

    @property
    def tsampling(self):
        return 1./self._frequency


class tmdataSR(tmdata):
    """tmdata with the memoization of rebinned data

    For a given cha and lbin, perform the rebin individually .
    """

    def __init__(self, filenames=[], channels=[], frequency=1):
        super(tmdataRebin, self).__init__(filenames, channels, frequency)
        self._rebindata = {}
        self._rebindata[1] = self.data
        self._rebinfreq = {}
        self._rebinfreq[1] = self.frequency

    def rebin(self, cha=1, lbin=1):
        self._checkchalbin(cha, lbin)
        try:
            if self._rebindata[lbin][cha-1] is not None:
                return self._rebindata[lbin][cha-1], self._rebinfreq[lbin]
            else:
                k, f = self._rebin(cha, lbin)
                self._rebindata[lbin][cha-1] = k
                return k, f
        except:
            self._rebindata.setdefault(lbin, [None]*self.nchannels)
            k, f = self._rebin(cha, lbin)
            self._rebindata[lbin][cha-1] = k
            self._rebinfreq[lbin] = f
            return k, f

    def _rebin(self, cha, lbin):
        for i in sorted(self._rebindata.keys(), reverse=True):
            if lbin % i == 0 and (self._rebindata[i][cha-1] is not None):
                temp_ds = self._rebindata[i][cha-1].size
                temp_lbin = lbin // i
                bds = temp_ds // (temp_lbin)
                k = (self._rebindata[i][cha-1][:bds*temp_lbin]
                            .reshape(bds, temp_lbin).sum(axis=1))
                f = self._rebinfreq[i] / temp_lbin
                return k, f

    @property
    def rebindata(self):
        return self._rebindata
    @property
    def rebinfrequency(self):
        return self._rebinfrequency


class tmdataMR(tmdata):
    """tmdata with the memoization of rebinned data

    For a given rebin-factor, this rebins data of all channels and memorizes
    them.
    """

    def __init__(self, filenames=[], channels=[], frequency=1):
        super(tmdataMR, self).__init__(filenames, channels, frequency)
        self._rebindata = {}
        self._rebindata[1] = self.data
        self._rebinfreq = {}
        self._rebinfreq[1] = self.frequency

    def rebin(self, cha=1, lbin=1):
        self._checkchalbin(cha, lbin)
        try:
            return self._rebindata[lbin][cha-1], self._rebinfreq[lbin]
        except:
            k, f = self._rebin(cha, lbin)
            self._rebindata[lbin] = k
            self._rebinfreq[lbin] = f
            return k[cha-1], f

    def _rebin(self, cha, lbin):
        for i in sorted(self._rebindata.keys(), reverse=True):
            if lbin % i == 0:
                temp_ds = self._rebindata[i][0].size
                temp_lbin = lbin // i
                bds = temp_ds // (temp_lbin)
                k = [temp_k[:bds*temp_lbin].reshape(bds, temp_lbin).sum(axis=1)
                        for temp_k in self._rebindata[i]]
                f = self._rebinfreq[i] / temp_lbin
                return k, f

    @property
    def rebindata(self):
        return self._rebindata
    @property
    def rebinfrequency(self):
        return self._rebinfrequency
