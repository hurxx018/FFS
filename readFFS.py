

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
            self.filenames = filenames
        else:
            raise TypeError("filenames is not a python list.")

        if isinstance(channels, list):
            self.channels = channels
        else:
            raise TypeError("channels is not a python list.")

        if isinstance(frequency, int):
            self.frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")


    def setFilenames(self, filenames):
        if isinstance(filenames, list):
            self.filenames = filenames
        else:
            raise TypeError("filenames is not a python list.")


    def setChannels(self, channels):
        if isinstance(channels, list):
            self.channels = channels
        else:
            raise TypeError("channels is not a python list.")


    def setFrequency(self, frequency):
        if isinstance(frequency, int):
            self.frequency = frequency
        else:
            raise TypeError("frequency is not an integer.")


    def getFilenames(self):
        return self.filenames


    def getChannels(self):
        return self.channels


    def getNChannels(self):
        return len(self.channels)


    def getFrequency(self):
        return self.frequency
