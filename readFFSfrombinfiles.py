import numpy as np
from readFFS import readFFS

class readFFSfrombinfiles(readFFS):
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
    def __init__(self, filenames=[], channels=[], frequency=0):
        readFFS.__init__(self, filenames, channels, frequency)
        self.data = self.readFFSData()

    def readFFSData(self):
        if self.getFilenames() == []:
            return []
        else:
            data = []
            for filename in self.getFilenames():
                with open(filename) as f:
                    temp_data = np.fromfile(f, dtype=np.int32)
                data.append(temp_data)
            return data


    def setFilenames(self, filenames):
        readFFS.setFilenames(self, filenames)
        self.data = self.readFFSData()


    def getData(self):
        return self.data


    def getInfo(self):
        """
        Returns the information in the class
        """
        return  {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}


    def __str__(self):
        for key, value in self.getInfo().items():
            print("{0}  :   {1}".format(key, value))
        return ""



def main():
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = readFFSfrombinfiles([filename1], [1])
    print("filenames : ")
    for x in ffsdata.getFilenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.getChannels())
    print("n_channels: ", ffsdata.getNChannels())
    print("data :")
    for data in ffsdata.getData():
        print(data[0:20])

    print("")

    ffsdata.setFilenames([filename1, filename2])
    ffsdata.setChannels([1, 2])
    print("filenames : ")
    for x in ffsdata.getFilenames():
        print("  {}".format(x))

    print("channels : ",ffsdata.getChannels())
    print("n_channels: ", ffsdata.getNChannels())
    print("data :")
    for data in ffsdata.getData():
        print(data[0:20])

    print("")

    ffsdata.setFilenames([filename1, filename2, filename3])
    ffsdata.setChannels([1, 2, 3])
    print("filenames : ")
    for x in ffsdata.getFilenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.getChannels())
    print("n_channels: ", ffsdata.getNChannels())
    print("data :")
    for data in ffsdata.getData():
        print(data[0:20])

if __name__=="__main__":
    main()
