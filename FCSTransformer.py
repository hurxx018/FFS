import itertools, collections
import numpy as np

from tmdata import tmdata
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
    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 4

    def __init__(self, segmentlength=DATACYCLE
                     , tau1_min = 1     # unrebinned data g(tau1_min : tau1_max)
                     , tau1_max = 16
                     , binfactor = 4    # rebin data consecutively by this factor
                     , taur_min = 4     # rebin data gr(taur_min : taur_max)
                     , taur_max = 16
                     , npts_limit = 64
                     , channels=[]):
        super(FCSTransformer, self).__init__(segmentlength, channels)
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

    # def fit(self, X, y=None):
    #     return self

    def transform(self, X):
        if not isinstance(X, tmdata):
            raise TypeError("The type of input is not tmdata.")
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
        time = self._time_index / (frequency * 1.)
        result[(0,0)] = time

        reshaped_data={}
        for x in itertools.combinations_with_replacement(channels, 2):
            arr1 = data[x[0] - 1]
            arr2 = data[x[1] - 1]
            result[x] = self.subCorrelations(arr1, arr2, x, reshaped_data)
        del reshaped_data  #delete reshaped_data
        return result

    def subCorrelations(self, arr1, arr2, indices, reshaped_data):
        n_segment = arr1.size // self._segmentlength
        x1, x2 = indices
        if x1 in reshaped_data:
            temp_arr1 = reshaped_data[x1]
        else:
            temp_arr1 = arr1[0:self._segmentlength*n_segment]  \
                            .reshape((n_segment, self._segmentlength))
            reshaped_data[x1] = temp_arr1

        if x2 in reshaped_data:
            temp_arr2 = reshaped_data[x2]
        else:
            temp_arr2 = arr2[0:self._segmentlength*n_segment]  \
                            .reshape((n_segment, self._segmentlength))
            reshaped_data[x2] = temp_arr2

        temp_shape = temp_arr1.shape

        correlations = np.zeros(len(self._time_index)) #[]
        correlations_stderr = np.zeros(len(self._time_index))#[]
        index = 0
        for t0 in range(self._tau1_min, self._tau1_max+1):
            t_arr1 = temp_arr1[:, 0:self._segmentlength-t0]
            t_arr2 = temp_arr2[:, t0:self._segmentlength]

            t_corr = (t_arr1*t_arr2).mean(axis=1)          \
                        / t_arr1.mean(axis=1)            \
                        / t_arr2.mean(axis=1) - 1.
            correlations[index] = t_corr.mean()
            correlations_stderr[index]   \
                    = t_corr.std() / np.sqrt(temp_shape[0] - 1)
            index += 1

        temp_segmentlength = self._segmentlength
        number_of_rebins = (np.log(self._segmentlength / self._npts_limit)  \
                        / np.log(self._binfactor) + 1).astype(int)
        for rbin in range(number_of_rebins - 1):
            temp_segmentlength = temp_segmentlength // self._binfactor
            y1, y2 = (x1, rbin), (x2, rbin)
            if y1 in reshaped_data:
                tt_arr1 = reshaped_data[y1]
            else:
                tt_arr1 = temp_arr1.reshape(temp_shape[0], temp_segmentlength, temp_arr1.shape[1]//temp_segmentlength).sum(axis=2)
                reshaped_data[y1] = tt_arr1
            if y2 in reshaped_data:
                tt_arr2 = reshaped_data[y2]
            else:
                tt_arr2 = temp_arr2.reshape(temp_shape[0], temp_segmentlength, temp_arr2.shape[1]//temp_segmentlength).sum(axis=2)
                reshaped_data[y2] = tt_arr2

            for t1 in range(self._taur_min, self._taur_max + 1):
                t_arr1 = tt_arr1[:, 0:temp_segmentlength - t1]
                t_arr2 = tt_arr2[:, t1:temp_segmentlength]

                t_corr = (t_arr1*t_arr2).mean(axis=1)           \
                                    / t_arr1.mean(axis=1)       \
                                    / t_arr2.mean(axis=1) - 1.
                correlations[index] = t_corr.mean()
                correlations_stderr[index]  = t_corr.std()      \
                                    / np.sqrt(temp_shape[0]-1)
                index += 1

        Correlations = collections.namedtuple("Correlations",      \
                            ["correlations", "correlations_stderr"])
        return Correlations(correlations, correlations_stderr)


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
    import matplotlib.pyplot as plt
    import time
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename1, filename2, filename3], \
                                        [1,2,3], frequency=100000)
    print("filenames : ")
    for x in ffsdata.filenames:
        print("  {}".format(x))

    print("channels : ", ffsdata.channels)

    cor = FCSTransformer(channels=[1, 2, 3])
    start = time.time()
    xx = cor.transform(ffsdata)
    end = time.time()
    print(end - start)
    plt.plot(xx[(0,0)], xx[(2,2)][0], label="2")
    plt.plot(xx[(0,0)], xx[(1,1)][0], label="1")
    plt.plot(xx[(0,0)], xx[(3,3)][0], label="3")
    plt.plot(xx[(0,0)], xx[(1,2)][0], label="12")
    plt.plot(xx[(0,0)], xx[(1,3)][0], label="13")
    plt.plot(xx[(0,0)], xx[(2,3)][0], label="23")
    plt.xscale("log")
    plt.legend()
    plt.show()
    # print(xx)

if __name__ == "__main__":
    main()
