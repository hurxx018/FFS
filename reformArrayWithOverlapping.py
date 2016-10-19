import numpy as np
#import
from readFFS import readFFS
from FFSTransformer import FFSTransformer



class reformArrayWithOverlapping(FFSTransformer):

    def __init__(self, segmentlength=32768, hoplength=None, channels=[], flatten=False):

        FFSTransformer.__init__(self, segmentlength, channels)
        if isinstance(hoplength, int):
            self._hoplength = hoplength
        elif isinstance(hoplength, float):
            self._hoplength = int(self.segmentlength*hoplength)
        elif hoplength == None:
            self._hoplength = self.segmentlength
        else:
            raise TypeError("The type of hoplength is either int or float.")
        self._flatten = flatten  # flatten an array of data for the other transformations such FCSTransformer.

    def transform(self, X):
        if self.channels == []:
            filenames = X.filenames[:]
            data = X.data
            channels = X.channels[:]
            frequency = X.frequency
            result = readFFS(filenames, channels, frequency)
            result.data = self.reformWithOverlapping(data, channels, frequency)
        else:
            filenames = [X.filenames[i-1] for i in self.channels]
            data = X.data
            channels = self.channels[:]
            frequency = X.frequency
            result = readFFS(filenames, channels, frequency)
            result.data = self.reformWithOverlapping(data, channels, frequency)
        return result

    def reformWithOverlapping(self, data, channels, frequency):
        sgl = self.segmentlength
        hl = self.hoplength
        result = []
        for i in channels:
            start = 0
            end = sgl
            temp = []
            while end <= data[i-1].size:
                try:
                    temp.append(data[i-1][start:end])
                    start += hl
                    end += hl
                except IndexError:
                    break
            if self._flatten:
                temp_forflattening = np.vstack(temp)
                result.append(np.flatten(temp_forflattening))
            else:
                result.append(np.vstack(temp))
        return result

    @property
    def hoplength(self):
        return self._hoplength

    @hoplength.setter
    def hoplength(self, value):
        self._hoplength = value

def main():
    from readFFSfrombinfiles import readFFSfrombinfiles
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    ffsdata = readFFSfrombinfiles([filename1, filename2], [1,2])
    reform_trans = reformArrayWithOverlapping(segmentlength=32768*1, hoplength=32768*4//4, channels=[2])
    temp = reform_trans.transform(ffsdata)
    print(temp.data[0].shape)
    print(temp.data[0].size)
    #print(temp)
    #print(temp.data[1].shape)
    #print(temp.data[1].size)
if __name__ == "__main__":
    main()
