import itertools
import numpy as np

from FFSTransformer import FFSTransformer

class PCH2DTransformer(FFSTransformer):
    """Calculate Photon Counting Histogram

    """
    def __init__(self, segmentlength=32768*4
                    , binfactor=[]
                    , channels=[]):
        FFSTransformer.__init__(self, segmentlength, channels)

        if binfactor == []:
            self._binfactor = [1]
        elif isinstance(binfactor, list):
            self._binfactor = binfactor
        elif isinstance(binfactor, int):
            self._binfactor = [binfactor]
        else:
            raise TypeError("The type of binfactor is not correct.")

    def transform(self, X):         #X = FFSdata
        if self._channels == []:
            data = X.data
            channels = X.channels
            frequency = X.frequency
            result = self.pch2destimators(data, channels, frequency)
        else:
            data = X.data
            frequency = X.frequency
            result = self.pch2destimators(data, self._channels, frequency)
        return result

    def pch2destimators(self, data, channels, frequency):
        result = {}
        for i in itertools.combinations(channels, 2):
            result[i] = self.sub_pch2destimators(data[i[0] - 1], data[i[1] - 1], frequency)

        return result

    def sub_pch2destimators(self, array1, array2, frequency):
        binf = self._binfactor[0]
        time = (self._segmentlength // binf)/float(frequency)

        n_segments = array1.size // self._segmentlength * binf
        temp_data1 = array1[0:n_segments*(self._segmentlength//binf)]
        # temp_data1 = (array1[0:n_segments*(self._segmentlength//binf)]
        #             .reshape((n_segments, self._segmentlength//binf)))
        temp_data2 = array2[0:n_segments*(self._segmentlength//binf)]
        # temp_data2 = (array2[0:n_segments*(self._segmentlength//binf)]
        #             .reshape((n_segments, self._segmentlength//binf)))

        kmin1, kmax1 = temp_data1.min(), temp_data1.max()
        kmin2, kmax2 = temp_data2.min(), temp_data2.max()
        kvec1 = np.arange(kmin1,kmax1+2)
        kvec2 = np.arange(kmin2,kmax2+2)
        pch, kvec = np.histogramdd([temp_data1, temp_data2], [kvec1, kvec2])

        temp_data1 = temp_data1.reshape((n_segments, self._segmentlength//binf))
        temp_data2 = temp_data2.reshape((n_segments, self._segmentlength//binf))

        # for loop is used << no idea for using apply_along_axis
        pchBySeg = np.array([np.histogramdd([temp_data1[i,:], temp_data2[i, :]], [kvec1, kvec2])[0] for i in range(n_segments)])

        # The current reture is temporary.
        return kvec, pch, pchBySeg


def main():
    import readFFSfrombinfiles as rffs
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename1, filename2], [1, 2], 100000)
    pch0 = PCH2DTransformer(segmentlength=32768*8)

    pch0_est = pch0.transform(ffsdata)

    print(pch0)
    print(pch0_est)
    print("---"*10)
if __name__ == "__main__":
    main()
