import numpy as np

from FFSTransformer import FFSTransformer

class PCH1DTransformer(FFSTransformer):
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
            result = self.pch1destimators(data, channels, frequency)
        else:
            data = X.data
            frequency = X.frequency
            result = self.pch1destimators(data, self._channels, frequency)
        return result

    def pch1destimators(self, data, channels, frequency):
        result = {}
        for i in channels:
            result[i] = self.sub_pch1destimators(data[i - 1], frequency)

        return result

    def sub_pch1destimators(self, array, frequency):
        binf = self._binfactor[0]
        time = (self._segmentlength // binf)/float(frequency)

        n_segments = array.size // self._segmentlength * binf
        temp_data = (array[0:n_segments*(self._segmentlength//binf)]
                    .reshape((n_segments, self._segmentlength//binf)))
        kmin = temp_data.min();
        kmax = temp_data.max();
        kvec = np.arange(kmin,kmax+2);
        pch = np.histogram(temp_data, bins=kvec)[0]
        nseg = temp_data.shape[0]
        #first attempt: using for loop
        #pchBySeg = np.zeros( (nseg, kmax+1), dtype=np.int64 )
        #for i in np.arange(25):
        #    pchBySeg[i,] = np.bincount(temp_data[i,],minlength=kmax+1)
        #pch = pchBySeg.sum(axis=0)
        #pchBySeg[:,0].std()
        #second approach using numpy
        """numpy.vstack is used to convert a list 1D arrays into 2D array instead of using two steps (tolist and array)."""
        pchBySeg = np.vstack(np.apply_along_axis(np.histogram,1,temp_data,bins=kvec)[:,0]) #np.array(np.apply_along_axis(np.histogram,1,temp_data,bins=kvec)[:,0].tolist())
        #pch = pchBySeg.sum(axis=0)

        # The current reture is temporary.
        return kvec, pch, pchBySeg


def main():
    import readFFSfrombinfiles as rffs
    filename2 = "A488_cal.2.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename2], [1], 100000)
    pch0 = PCH1DTransformer(segmentlength=32768*8)

    pch0_est = pch0.transform(ffsdata)

    print(pch0)
    print(pch0_est)
    print("---"*10)
if __name__ == "__main__":
    main()
