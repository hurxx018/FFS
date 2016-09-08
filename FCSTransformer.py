import itertools, collections
import numpy as np

from FFSTransformer import FFSTransformer

class FCSTransformer(FFSTransformer):
    """
    Calculate auto-correlations and cross-correlations of photon count data

    Default:
    segmentlength : 32768L*4
    channels : [] <<== consider all channels in the data,
        otherwise some channels can be selected,
        for example,
        [1] for the first channel
        [1, 2] for the first and second channels and so on
    """
    def __init__(self, segmentlength=32768*4
                     , tau1_min = 1     # unrebinned data g(tau1_min : tau1_max)
                     , tau1_max = 16
                     , binfactor = 4    # rebin data consecutively by this factor
                     , taur_min = 4     # rebin data gr(taur_min : taur_max)
                     , taur_max = 16
                     , npts_limit = 64
                     , channels=[]):

        FFSTransformer.__init__(self, segmentlength, channels)
        self._tau1_min = tau1_min
        self._tau1_max = tau1_max
        self._binfactor = binfactor
        self._taur_min = taur_min
        self._taur_max = taur_max
        self._npts_limit = npts_limit
        self._time_index = self.generate_time_index()

    def generate_time_index(self):
        time_index = [np.arange(self.tau1_min, self.tau1_max + 1)]
        number_of_rebins = (np.log(self.segmentlength/self.npts_limit)      \
                                / np.log(self.binfactor) + 1).astype(int)
        for i in range(1, number_of_rebins):
            time_index.append(np.arange(self.taur_min, self.taur_max + 1)   \
                                        * self.binfactor**i)
        return np.concatenate(time_index)

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        if self.channels == []:
            data = X.data
            channels = X.channels
            frequency = X.frequency
            result = self.calCorrelations(data, channels, frequency)
        else:
            data = X.data
            frequency = X.frequency
            result = self.calCorrelations(data, self.channels, frequency)
        return result

    def calCorrelations(self, data, channels, frequency):
        result = {}
        for x in itertools.combinations_with_replacement(channels, 2):
            arr1 = data[x[0] - 1]
            arr2 = data[x[1] - 1]
            result[x] = self.subCorrelations(arr1, arr2, frequency)
        return result

    def subCorrelations(self, arr1, arr2, frequency):
        assert arr1.size == arr2.size

        time = self._time_index / (frequency * 1.)
        n_segment = arr1.size // self.segmentlength
        temp_arr1 = arr1[0:self.segmentlength*n_segment]  \
                            .reshape((n_segment, self.segmentlength))
        temp_arr2 = arr2[0:self.segmentlength*n_segment]  \
                            .reshape((n_segment, self.segmentlength))
        temp_shape = temp_arr1.shape

        correlations = np.zeros(len(time)) #[]
        correlations_stderr = np.zeros(len(time))#[]
        index = 0
        for t0 in range(self.tau1_min, self.tau1_max+1):
            t_arr1 = temp_arr1[:, 0:self.segmentlength-t0]*1.
            t_arr2 = temp_arr2[:, t0:self.segmentlength]*1.

            t_corr = np.mean(t_arr1*t_arr2, axis=1)          \
                        / np.mean(t_arr1, axis=1)            \
                        / np.mean(t_arr2, axis=1) - 1.
            correlations[index] = np.mean(t_corr)
            correlations_stderr[index]   \
                    = np.std(t_corr) / np.sqrt(temp_shape[0] - 1)
            index += 1

        temp_segmentlength = self.segmentlength
        number_of_rebins = (np.log(self.segmentlength / self.npts_limit)  \
                        / np.log(self.binfactor) + 1).astype(int)
        for rbin in range(number_of_rebins - 1):
            temp_segmentlength = temp_segmentlength // self.binfactor
            tt_arr1 = self.rebin(temp_arr1 * 1., \
                                    (temp_shape[0], temp_segmentlength))
            tt_arr2 = self.rebin(temp_arr2 * 1., \
                                    (temp_shape[0], temp_segmentlength))

            for t1 in range(self.taur_min, self.taur_max + 1):
                t_arr1 = tt_arr1[:, 0:temp_segmentlength - t1]
                t_arr2 = tt_arr2[:, t1:temp_segmentlength]

                t_corr = np.mean(t_arr1*t_arr2, axis=1)             \
                                    / np.mean(t_arr1, axis=1)       \
                                    / np.mean(t_arr2, axis=1) - 1.
                correlations[index] = np.mean(t_corr)
                correlations_stderr[index]  = np.std(t_corr)           \
                                    / np.sqrt(temp_shape[0]-1)
                index += 1

        Correlations = collections.namedtuple("Correlations",      \
                            ["time","correlations", "correlations_stderr"])
        return Correlations(time, correlations, correlations_stderr)


    def rebin(self, a, newshape):
        """
        Rebin a FFS data
        """
        M, N = a.shape
        m, n = newshape
        if n<N:
            t0 = a.reshape((m, M // m, n, N // n))
            return np.sum(np.sum(t0, axis=3), axis=1)
        else:
            return np.repeat(np.repeat(a, m // M, axis=0), n // N, axis=1)

    @property
    def tau1_min(self):
        return self._tau1_min

    @tau1_min.setter
    def tau1_min(self, value):
        self._tau1_min = value
        self._time_index = self.generate_time_index()

    @property
    def tau1_max(self):
        return self._tau1_max

    @tau1_max.setter
    def tau1_max(self, value):
        self._tau1_max = value
        self._time_index = self.generate_time_index()

    @property
    def binfactor(self):
        return self._binfactor

    @binfactor.setter
    def binfactor(self, value):
        self._binfactor = value
        self._time_index = self.generate_time_index()

    @property
    def taur_min(self):
        return self._taur_min

    @taur_min.setter
    def taur_min(self, value):
        self._taur_min = value
        self._time_index = self.generate_time_index()

    @property
    def taur_max(self):
        return self._taur_max

    @taur_max.setter
    def taur_max(self, value):
        self._taur_max = value
        self._time_index = self.generate_time_index()

    @property
    def npts_limit(self):
        return self._npts_limit

    @npts_limit.setter
    def npts_limit(self, value):
        self._npts_limit = value


def main():
    import readFFSfrombinfiles as rffs
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename1, filename2, filename3], \
                                        [1,2,3], frequency=100000)
    print("filenames : ")
    for x in ffsdata.filenames:
        print("  {}".format(x))

    print("channels : ", ffsdata.channels)

    cor = FCSTransformer(channels=[1, 2])
    print(cor)
    xx = cor.transform(ffsdata)
    print(xx)

if __name__ == "__main__":
    main()
