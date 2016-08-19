import itertools, collections
import numpy as np

class calfcsTransformer():
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
        self._segmentlength = segmentlength
        self._tau1_min = tau1_min
        self._tau1_max = tau1_max
        self._binfactor = binfactor
        self._taur_min = taur_min
        self._taur_max = taur_max
        self._npts_limit = npts_limit
        self._channels = channels

    def fit(self, X, y=None):
        return self

    def transform(self, X):

        if self._channels == []:
            data = X.get_data()
            channels = X.get_channels()
            result = self.cal_correlations(data, channels)
        else:
            data = X.get_data()
            result = self.cal_correlations(data, self._channels)
        return result

    def cal_correlations(self, data, channels):
        result = {}
        for x in itertools.combinations_with_replacement(channels, 2):
            arr1 = data[x[0]-1]
            arr2 = data[x[1]-1]
            result[x] = self.sub_correlations(arr1, arr2)
        return result

    def sub_correlations(self, arr1, arr2):
        assert arr1.size == arr2.size

        n_segment = arr1.size//self._segmentlength
        temp_arr1 = arr1[0:self._segmentlength*n_segment].reshape((n_segment, self._segmentlength))
        temp_arr2 = arr2[0:self._segmentlength*n_segment].reshape((n_segment, self._segmentlength))
        temp_shape = temp_arr1.shape

        correlations = []
        correlations_stderr = []
        for t0 in range(self._tau1_min, self._tau1_max+1):
            t_arr1 = temp_arr1[:, 0:self._segmentlength-t0]*1.
            t_arr2 = temp_arr2[:, t0:self._segmentlength]*1.

            t_corr = np.mean(t_arr1*t_arr2, axis=1)/np.mean(t_arr1, axis=1)/np.mean(t_arr2, axis=1)-1.
            correlations.append(np.mean(t_corr))
            correlations_stderr.append(np.std(t_corr)/np.sqrt(temp_shape[0]-1))

        temp_segmentlength = self._segmentlength
        number_of_rebins = int(np.log(self._segmentlength/self._npts_limit)/np.log(self._binfactor) + 1)
        for rbin in range(number_of_rebins):
            temp_segmentlength = temp_segmentlength//self._binfactor
            tt_arr1 = self.rebin(temp_arr1*1., (temp_shape[0], temp_segmentlength))
            tt_arr2 = self.rebin(temp_arr2*1., (temp_shape[0], temp_segmentlength))

            for t1 in range(self._taur_min, self._taur_max+1):
                t_arr1 = tt_arr1[:, 0:temp_segmentlength-t1]
                t_arr2 = tt_arr2[:, t1:temp_segmentlength]

                t_corr = np.mean(t_arr1*t_arr2, axis=1)/np.mean(t_arr1, axis=1)/np.mean(t_arr2, axis=1) - 1.
                correlations.append(np.mean(t_corr))
                correlations_stderr.append(np.std(t_corr)/np.sqrt(temp_shape[0]-1))

        Correlations = collections.namedtuple("Correlations", ["correlations", "correlations_stderr"])
        return Correlations(correlations, correlations_stderr)

    def rebin(self, a, newshape ):
        M, N = a.shape
        m, n = newshape
        if n<N:
            t0 = a.reshape((m,M//m,n,N//n))
            return np.mean(np.mean(t0, axis=3), axis=1)
        else:
            return np.repeat(np.repeat(a, m//M, axis=0), n//N, axis=1)

    def get_params(self):
        return {key:value for key, value in self.__dict__.items()
                         if not key.startswith('__') and not callable(key)}


def main():
    import readFFSfrombinfiles as rffs
    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename1, filename2, filename3], [1,2,3])
    print("filenames : ")
    for x in ffsdata.get_filenames():
        print("  {}".format(x))

    print("channels : ", ffsdata.get_channels())

    cor = calfcsTransformer(channels=[1,2])
    print(cor.get_params())
    # xx = cor.transform(ffsdata)
    # print(xx.get_parmas())

if __name__ == "__main__":
    main()
