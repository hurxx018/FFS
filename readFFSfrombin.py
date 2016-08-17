import os
import numpy as np

class readFFSfrombinfiles(object):
    """
    Read FFS data from bin files (###.bin) which are collected at timemode
    input variables:
    filenames : a python list of filenames
    channels : the number of channels in FFS data.

    output variables;
    data : a python list of numpy arrays.

    instance methods:
    readFFSData : read a data for each filename (32bit integer)
    """
    def __init__(self, filenames=[], channels=[]):
        if isinstance(filenames, list):
            self._filenames = filenames
        else:
            raise TypeError("filenmes is not a python list.")

        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")

        self._data = self.readFFSData()

    def readFFSData(self):
        if self._filenames == []:
            return []
        else:
            data = []
            for filename in self._filenames:
                with open(filename) as f:
                    temp_data = np.fromfile(f, dtype=np.int32)
                data.append(temp_data)
            return data

    def set_filenames(self, filenames):
        if isinstance(filenames, list):
            self._filenames = filenames
        else:
            raise TypeError("filenmes is not a python list.")

        self._data = self.readFFSData()

    def set_channels(self, channels):
        if isinstance(channels, list):
            self._channels = channels
        else:
            raise TypeError("channels is not a python list.")

    def get_filenames(self):
        return self._filenames

    def get_channels(self):
        return self._channels

    def get_nchannels(self):
        return len(self._channels)

    def get_data(self):
        return self._data


def main():
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = readFFSfrombinfiles([filename1], [1])
    print("filenames : ")
    for x in ffsdata.get_filenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.get_channels())
    print("n_channels: ", ffsdata.get_nchannels())
    print("data :")
    for data in ffsdata.get_data():
        print(data[0:20])

    print("")

    ffsdata.set_filenames([filename1, filename2])
    ffsdata.set_channels([1, 2])
    print("filenames : ")
    for x in ffsdata.get_filenames():
        print("  {}".format(x))

    print("channels : ",ffsdata.get_channels())
    print("n_channels: ", ffsdata.get_nchannels())
    print("data :")
    for data in ffsdata.get_data():
        print(data[0:20])

    print("")

    ffsdata.set_filenames([filename1, filename2, filename3])
    ffsdata.set_channels([1, 2, 3])
    print("filenames : ")
    for x in ffsdata.get_filenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.get_channels())
    print("n_channels: ", ffsdata.get_nchannels())
    print("data :")
    for data in ffsdata.get_data():
        print(data[0:20])

if __name__=="__main__":
    main()
