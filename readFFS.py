

class readFFS(object):
    """
    A superclass for the FFS data
    input variables:
    filenames : a python list of filenames
    channels : a python list of channels in FFS data.
    frequency: an integer of frequency of FFS data collection

    """

    def __init__(self, filenames=[], channels=[], frequency=0):
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


    def setFilenames(self, filenames):
        if isinstance(filenames, list):
            self._filenames = filenames
        else:
            raise TypeError("filenames is not a python list.")


    def setChannels(self, channels):
        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")


    def setFrequency(self, frequency):
        if isinstance(frequency, int):
            self._frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")


    def getFilenames(self):
        return self._filenames


    def getChannels(self):
        return self._channels


    def getNChannels(self):
        return len(self._channels)


    def getFrequency(self):
        return self._frequency
