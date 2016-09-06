

class readFFS(object):
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
    def nchannels(self):
        return len(self._channels)
