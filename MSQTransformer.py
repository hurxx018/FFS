import collections
import numpy as np

from FFSTransformer import FFSTransformer

class MSQTransformer(FFSTransformer):
    """Calculate mean of segmented Q values

    """
    def __init__(self, segmentlength=32768*4
                     , binfactors=[]
                     , channels=[]
                     , tshift=0):
        FFSTransformer.__init__(self, segmentlength, channels)

        if binfactors == []:
            self._binfactors = 2**np.arange(14)
        else:
            self._binfactors = np.array(binfactors).astype(int)
        self._tshift = tshift

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        if self._channels == []:
            data = X.data
            channels = X.channels
            frequency = X.frequency
            result = self.calqestimators(data, channels, frequency)
        else:
            data = X.data
            channels = X.channels
            frequency = X.frequency
            result = self.calqestimators(data, self._channels, frequency)
        return result

    def calqestimators(self, data, channels, frequency):
        result = {}
        for i in channels:
            if self._tshift == 0:
                result[i] = self.sub_qestimators(data[i - 1], frequency)
            else:
                result[i] = self.sub_qestimators_withtshift(data[i - 1], \
                                                    frequency)
        return result

    def sub_qestimators(self, array, frequency):
        time = np.zeros(self._binfactors.size)
        qestimator = np.zeros(self._binfactors.size)
        qestimator_stderr = np.zeros(self._binfactors.size)
        for i, binf in enumerate(self._binfactors):
            time[i] = (self._segmentlength // binf)/float(frequency)

            n_segments = array.size // self._segmentlength * binf
            temp_data = array[0:n_segments*(self._segmentlength//binf)] \
                        .reshape((n_segments, self._segmentlength//binf))
            temp_k1 = np.mean(temp_data, axis=1)
            temp_k2 = np.mean((temp_data*temp_data), axis=1)
            index = temp_k1 > 0.
            temp_q = temp_k2[index]/temp_k1[index] - temp_k1[index] - 1.
            qestimator[i] = np.mean(temp_q)
            if n_segments > 1:
                qestimator_stderr[i] = np.std(temp_q)/np.sqrt(n_segments - 1)
            else:
                qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1 #TO DO

        qestimator_tuple = collections.namedtuple("Qestimator", \
                    ["time", "qestimator", "qestimator_stderr", "tsfhit"])
        return qestimator_tuple(time, qestimator, qestimator_stderr, self._tshift)

    def sub_qestimators_withtshift(self, array, frequency):
        time = np.zeros(self._binfactors.size)
        qestimator = np.zeros(self._binfactors.size)
        qestimator_stderr = np.zeros(self._binfactors.size)
        index = 0
        for i, binf in enumerate(self._binfactors):
            time[i] = (self._segmentlength // binf)/float(frequency)

            n_segments = array.size // self._segmentlength * binf
            temp_data = (array[0:n_segments*(self._segmentlength//binf)]
                        .reshape((n_segments, self._segmentlength//binf)))

            mean1 = (np.mean(temp_data[:, :-self._tshift], axis=1)
                       .reshape((-1, 1)))
            mean2 = (np.mean(temp_data[:, self._tshift:], axis=1)
                       .reshape((-1, 1)))
            temp_q = (np.mean((temp_data[:, :-self._tshift] - mean1)
                        *(temp_data[:, self._tshift:] - mean2), axis=1)
                        .reshape((-1, 1))
                        /(mean1 + mean2) * 2 )
            qestimator[i] = np.mean(temp_q)
            if n_segments > 1:
                qestimator_stderr[i] = np.std(temp_q)/np.sqrt(n_segments - 1)
            else:
                qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1

        qestimator_tuple = collections.namedtuple("Qestimator", \
                    ["time", "qestimator", "qestimator_stderr", "tshift"])
        return qestimator_tuple(time, qestimator, qestimator_stderr, self._tshift)


    @property
    def binfactors(self):
        return self._binfactors

    @binfactors.setter
    def binfactors(self, values):
        if isinstance(values, np.ndarray):
            self._binfactors = values
        else:
            self._binfactors = np.array(values)

    @property
    def tshift(self):
        return self._tshift

    @tshift.setter
    def tshift(self, value):
        self._tshift = value


def main():
    import readFFSfrombinfiles as rffs
    filename2 = "A488_cal.2.001.bin"
    ffsdata = rffs.readFFSfrombinfiles([filename2], [1], 100000)
    msq0 = MSQTransformer(segmentlength=32768*8, tshift=0)
    msq1 = MSQTransformer(segmentlength=32768*8, tshift=1)

    msq0_est = msq0.transform(ffsdata)
    msq1_est = msq1.transform(ffsdata)

    print(msq0)
    print(msq0_est)
    print("---"*10)
    print(msq1)
    print(msq1_est)
if __name__ == "__main__":
    main()
