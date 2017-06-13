import numpy as np
from tmdata import tmdataSR

def iss0(filenames, nch):
    data = []
    for filename in self.filenames:
        with open(filename) as f:
            temp_data = np.fromfile(f, dtype=np.int32)
        data.append(temp_data)
    return data

def flex8(filenames, nch):
    with open(self.filenames[0]) as f:
        temp_data = np.fromfile(f, dtype=np.uint8)
    if nch == 1:
        return [temp_data]
    elif nch == 2:
        temp_data = (temp_data.reshape(temp_data.size//2, 2)
                                      .T)
        return [temp_data[0], temp_data[1]]
    else:
        raise ValueError("channels are incorrectly assigned.")

def flex16(filenames, nch):
    with open(self.filenames[0]) as f:
        temp_data = np.fromfile(f, dtype=np.uint16)
    if nch == 1:
        return [temp_data]
    elif nch == 2:
        temp_data = (temp_data.reshape(temp_data.size//2, 2)
                                      .T)
        return [temp_data[0], temp_data[1]]
    else:
        raise ValueError("channels are incorrectly assigned.")


class rFFS(tmdataSR):
    """
    Read FFS data from bin files (###.bin) which are collected at timemode
    input variables:
    filenames : a python list of filenames
    channels : a python list of channels in FFS data.
    frequency: an integer of frequency of FFS data collection

    output variables;
    data : a python list of numpy arrays.

    instance methods:
    readFFSData : read a data for each filename (32bit integer)
    """
    daqcards = {"iss0":iss0, "flex8":flex8, "flex16":flex16}

    def __init__(self, filenames=[], channels=[], frequency=1, daqcard=None):
        super(rFFS, self).__init__(filenames, channels, frequency)
        #tmdata.__init__(self, filenames, channels, frequency)
        self._data = self.readFFSData()
        self._checkdata()

    def readFFSData(self):
        if not self.filenames:
            return []
        else:
            return self.daqcards[daqcard](self.filenames, self.nchannels)

    @property
    def filenames(self):
        return self._filenames

    @filenames.setter
    def filenames(self, filenames):
        if isinstance(filenames, list):
            self._filenames = filenames
            self._data = self.readFFSData()
            self._checkdata()
        else:
            raise TypeError("filenames is not a python list.")

    def info(self):
        """
        Returns the information in the class
        """
        return  {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}


    def __str__(self):
        for key, value in self.info().items():
            print("{0}  :   {1}".format(key, value))
        return ""
