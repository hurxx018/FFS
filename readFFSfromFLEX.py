import numpy as np
from tmdata import tmdata


class readFFSfromFLEX(tmdata):
    """Read FFS data from dat files obtained by the FLEX card.

    input variables:
    filename : a list of a single filename
    channels : a python list of channels in FFS data.
    ex) [1] for single-channel FFS, [1, 2] for dual-channel FFS
    frequency: an integer of frequency of FFS data collection
    bit : 8 for 8 bit mode (Default), 16 for 16 bit mode

    output variables;
    data : a python list of numpy arrays.





    """
    def __init__(self, filenames=[], channels=[], frequency=1, bit=8):
        super(readFFSfromFLEX, self).__init__(filenames, channels, frequency)
        #readFFS.__init__(self, filenames, channels, frequency)
        self._bit = bit
        self._data = self.readFFSData()
        self._checkdata()

    def readFFSData(self):
        if self.filenames == []:
            return []
        else:
            if self._bit == 8:
                datatype = np.uint8
            elif self._bit == 16:
                datatype = np.uint16
            else:
                raise ValueError("bit is not correct >> bit is either 8 or 16")

            with open(self.filenames[0]) as f:
                temp_data = np.fromfile(f, dtype=datatype)

            if len(self.channels) == 1:
                data = [temp_data]
            elif len(self.channels) == 2:
                temp_data = (temp_data.reshape(temp_data.size//2, 2)
                                      .T)
                data = [temp_data[0], temp_data[1]]
            else:
                raise ValueError("channels are incorrectly assigned.")
            return data

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

    @property
    def bit(self):
        return self._bit
    @bit.setter
    def bit(self, bit_value):
        if bit_value in [8, 16]:
            self._bit = bit_value
        else:
            raise ValueError("bit is not correct. >> bit is either 8 or 16")

    def info(self):
        """
        Returns the information in the class
        """
        return  {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}

    def __str__(self):
        for key, value in self.info().items():
            print("{0}  :   {1}".format(key[1:], value))
        return ""


def main():
    filename1 = "alexa_05.dat"
    ffsdata = readFFSfromFLEX([filename1], [1, 2], 100000, 8)
    print(ffsdata)

    print("data :")
    for data in ffsdata.data:
        print(data[0:20])

if __name__ == "__main__":
    main()
