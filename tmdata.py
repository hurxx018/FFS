import numpy as np

class tmdata(object):
    """
    A superclass for the FFS data   v 0.0.1 in Python 3
    input variables:
    channels : a python list of channels in FFS data.
    frequency: an integer of frequency of FFS data collection
    memorize : a flag varialbe which enables to memorize rebinned data

    Instance Methods:
    Setters: setFilenames, setChannels, setFrequency
    Getters: getFilenames, getChannels, getNChannels, getFrequency

    """

    def __init__(self, channels=[], frequency=1, memorize=True):

        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")

        if isinstance(frequency, int):
            self._frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")
        self._data = {}
        self._filenames = {}
        self.memorize = memorize

    def rebin(self, ch=1, lbin=1):
        self._checkchalbin(ch, lbin)
        k = self._data[ch].get(lbin, None)
        if k is None:
            k = self._newrebin(ch, lbin)
            if self.memorize:
                self._data[ch][lbin] = k
        return k

    def _newrebin(self, ch, lbin):
        for i in sorted(self._data[ch].keys(), reverse=True):
            if lbin % i == 0:
                temp_ds = self._data[ch][i].size
                temp_lbin = lbin // i
                bds = temp_ds // (temp_lbin)
                k = (self._data[ch][i][:bds*temp_lbin]
                            .reshape(bds, temp_lbin).sum(axis=1))
                return k

    def _checkdata(self):
        if self.nchannels > 1:
            length = self._data[self._channels[0]][1].size
            for key, data in self._data.items():
                assert data[1].size == length, \
                    "The number of data points does not match."

    def _checkchalbin(self, ch, lbin):
        if ch not in self._channels:
            raise ValueError("{} is not available".format(ch))
        if not isinstance(lbin, int):
            raise TypeError('The type of lbin is int.')

    def info(self):
        """
        Returns the information in the class
        """
        return  {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def reset(self):
        if self._data:
            self._data = { key: {1: d[1]} for key, d in self._data.items()}

    def __str__(self):
        for key, value in self.info().items():
            print("{0}  :   {1}".format(key, value))
        return ""

    def ch(self, ch=None, lbin=1):
        if not ch:
            for i in sorted(self._channels):
                return self.rebin(ch=i, lbin=lbin)
        return self.rebin(ch=ch, lbin=lbin)

    @property
    def filenames(self):
        return self._filenames

    @property
    def channels(self):
        return self._channels

    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, frequency):
        # reset the self._frequency and assign a new frequency
        if isinstance(frequency, (int, float)):
            self._frequency = frequency
        else:
            raise TypeError("frequency is not a numeric.")

    # get a frequency after rebin with lbin
    def getfreq(self, lbin=1):
        return self._frequency*1./lbin

    @property
    def data(self):
        # get the original data only.
        return { key: d[1] for key, d in self._data.items() }

    @data.setter
    def data(self, data):
        if not isinstance(data, list):
            raise TypeError("The type of data is list.")
        if len(data) != self.nchannels:
            print("The number of channels is required to be {}".format(self.nchannels))
            raise ValueError("The number of channels is not correct")
        self._data = { i: {1:d} for i, d in zip(self._channels, data) }

    # not sure yet the way to deal with the filename because
    # flex card needs only one filename but iss or pp-1000 card needs multiple filenames
    # for multiple channels
    def setdata(self, ch, data, filename=None):
        if not isinstance(data, (np.ndarray)):
            raise TypeError("The type of data has tobe a numpy.ndarray.")
        if ch not in self._channels:
            print("The channel has to be in {}".format(self._channels))
            raise ValueError("The channel is not correct")
        self._data[ch] = {1: data}
        self._filenames[ch] = filename

    @property
    def nchannels(self):
        return len(self._channels)

    @property
    def tsampling(self):
        return 1./self._frequency[1]

def main():
    import numpy as np

    temp = tmdata([0], 20000)
    temp.setdata(0, np.ones(32768*2))
    temp.ch()


    temp = tmdata(['A', 'B'], 100000)
    temp.setdata('A', np.ones(32768*2))
    temp.setdata('B', np.ones(32768*2))
    print(temp)
    temp.ch('A',lbin=5)


    from rFFS import rFFS
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    rffs = rFFS( [1, 2, 3], 100000, "iss0" )
    ffsdata = rffs.load( [filename1, filename2, filename3] )
    print(ffsdata)
    print("data :")
    for key, data in ffsdata.data.items():
       print(key, ": ", data[0:20])
    print("")

    print("Rebin by factor 4")
    rebinned_k1 = ffsdata.rebin(1, 4)  # for the channel 1

    print("After rebinning, the _rebindata has a key of 4 and its rebinned data.")
    print(ffsdata)
    print()

    print("Rebin by factor 8")
    rebinned_k18 = ffsdata.rebin(1, 8)  # for the channel 1

    print("Rebin by factor 8")
    rebinned_k182 = ffsdata.rebin(1, 8)  # for the channel 1

    rebinned_k2 = ffsdata.rebin(2, 4)  # for the channel 2
    print("After rebinning, the _rebindata has a key of 4 and its rebinned data.")
    print(ffsdata)
    print()

    # rebin channel 3 with a rebin-factor 16
    rebinned_k3 = ffsdata.rebin(3, 16)  # for the channel 3
    print("After rebinning, the _rebindata has a key of 16 and its rebinned data.")
    print(ffsdata)
    print()


if __name__ == '__main__':
    main()
