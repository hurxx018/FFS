import numpy as np
from tmdata import tmdata, tmdataSR, tmdataMR

def iss0(filenames, nch):
    data = []
    for filename in filenames:
        with open(filename) as f:
            temp_data = np.fromfile(f, dtype=np.int32)
        data.append(temp_data)
    return data

def flex8(filenames, nch):
    with open(filenames[0]) as f:
        temp_data = np.fromfile(f, dtype=np.uint8)
    if nch == 1:
        return [temp_data]
    elif nch == 2:
        temp_data = (temp_data.reshape(temp_data.size//2, 2).T)
        return [temp_data[0], temp_data[1]]
    else:
        raise ValueError("channels are incorrectly assigned.")

def flex16(filenames, nch):
    with open(filenames[0]) as f:
        temp_data = np.fromfile(f, dtype=np.uint16)
    if nch == 1:
        return [temp_data]
    elif nch == 2:
        temp_data = (temp_data.reshape(temp_data.size//2, 2).T)
        return [temp_data[0], temp_data[1]]
    else:
        raise ValueError("channels are incorrectly assigned.")


class rFFS(object):
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
    daqcards = {"iss0":iss0, "flex8":flex8, "flex16":flex16, "sim":0}

    def __init__(self, channels=[], frequency=1, daqcard=None):
        if daqcard not in self.daqcards.keys():
            print("daqcard can be either one among")
            for key in self.daqcards.keys():
                print(key)
            raise ValueError("daqcard is not given correctly.")
        self._daqcard = daqcard
        self._channels = channels
        self._frequency = frequency

    def load(self, filenames=[]):
        tm_dataSR = tmdataSR(filenames, self._channels, self._frequency)
        if filenames:
            if self._daqcard != "sim":
                tm_dataSR.data = self.daqcards[self._daqcard](filenames, tm_dataSR.nchannels)
                tm_dataSR._checkdata()
        tm_dataSR.assign_rebindata()
        return tm_dataSR

    # def get_tmdataMR(self, tm_data):     This can be deprecatted for the easy management
    #     # get a tmdataMR for the convenience of rebinning
    #     if not isinstance(tm_data, tmdata):
    #         raise TypeError("Type of input is not tmdata.")
    #     tm_dataMR = tmdataMR(tm_data.filenames, tm_data.channels, tm_data.frequency)
    #     tm_dataMR.data = tm_data.data
    #     tm_dataMR.assign_rebindata()
    #     return tm_dataMR

def main():

    # Use the iss cards
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = rFFS( [1, 2, 3], 100000, "iss0" )
    temp = ffsdata.load( [filename1, filename2, filename3] )
    print(temp)
    print("data :")
    for data in temp.data:
        print(data[0:20])
    print("")

    # Use the Flex card
    filename1 = "alexa_05.dat"
    ffsdata = rFFS([1, 2], 100000, "flex8") # 8bit
    temp = ffsdata.load([filename1])
    print(temp)
    print("data :")
    for data in temp.data:
        print(data[0:20])
    print()

    # For the simulation
    ffsdata = rFFS([1, 2], 100000, "sim")
    temp = ffsdata.load([])
    temp.data = [np.ones(32768*2, dtype=int), np.ones(32768*2, dtype=int)]
    temp.assign_rebindata()
    print(temp)
    for data in temp.data:
        print(data[0:20])

if __name__ == '__main__':
    main()
