import collections
import numpy as np

from tmdata import tmdata
from FFSTransformer import FFSTransformer

def dark(X_dark, binfactors, segmentlength=32768*4, channel=1, timeshift=1):
    """ Calculate k1 and k2 from tmdata in terms of binfactors

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    timeshift >>  a factor of timeshift; an integer

    ...

    """
    k1 = np.zeros( len(binfactors) )
    k2 = np.zeros( len(binfactors) )
    data = X_dark.ch(channel, 1)
    for i, b in enumerate(binfactors):
        nsegs = data.size//segmentlength*b
        temp = data[:nsegs*segmentlength//b].reshape(nsegs, segmentlength//b)
        if timeshift > 0:
            t_k1_1 = data[:, :-timeshift].mean(axis=1)
            t_k1_2 = data[:, timeshift:].mean(axis=1)
            t_k2 = (data[:, :-timeshift]*data[:, timeshift:]).mean(axis=1) - t_k1_1*t_k1_2
            k1[i] = (t_k1_1 + t_k1_2).mean(axis=0)/2.
            k2[i] = t_k2.mean(axis=0)
        elif timeshift == 0:
            t_k1 = data.mean(axis=1)
            t_k2 = np.power(data, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)
            k1[i] = t_k1.mean(axis=0)
            k2[i] = t_k2.mean(axis=0)
        else:
            raise ValueError("timeshift is an integer of 0 or positive.")
    return (k1, k2)


def msq_est_dark(X, binfactors, segmentlength=32768*4, channel=1, timeshift=1, dark=None):
    """ Calculate MSQ estimators with dark count corrections

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    timeshift >>  a factor of timeshift; an integer

    dark >> [k1, k2] from a function dark

    """
    if dark:
        assert len(binfactors) == dark[0].size, "inconsistency in the size of dark k1 and k2"
    else:
        raise ValueError("dark is not available.")

    msq = np.zeros( len(binfactors) )
    msq_err = np.zeros( len(binfactors) )

    data = X.ch(channel, 1)
    for i, b in enumerate(binfactors):
        nsegs = data.size//segmentlength*b

        dark1 = np.ones(nsegs) * dark[0][i]
        dark2 = np.ones(nsegs) * dark[1][i]

        temp = data[:nsegs*segmentlength//b].reshape(nsegs, segmentlength//b)
        if timeshift > 0:
            t_k1_1 = data[:, :-timeshift].mean(axis=1)
            t_k1_2 = data[:, timeshift:].mean(axis=1)
            t_k2 = (data[:, :-timeshift]*data[:, timeshift:]).mean(axis=1) - t_k1_1*t_k1_2

            m1 = (t_k1_1 + t_k1_2)/2. - dark1   # dark count correction
            m2 = t_k2 - dark2                   # dark count correction
            q = m2/m1
            msq[i] = q.mean()
            msq_err[i] = q.std()

        elif timeshift == 0:
            t_k1 = data.mean(axis=1)
            t_k2 = np.power(data, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)

            m1 = t_k1 - dark1
            m2 = t_k2 - dark2
            q = m2/m1
            msq[i] = q.mean()
            msq_err[i] = q.std()

        else:
            raise ValueError("timeshift is an integer of 0 or positive.")
    return (msq, msq_err)


def msq_est(X, binfactors, segmentlength=32768*4, channel=1, timeshift=1):
    """Calculate MSQ estimators without dark count corrections

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    timeshift >>  a factor of timeshift; an integer

    """

    msq = np.zeros( len(binfactors) )
    msq_err = np.zeros( len(binfactors) )

    data = X.ch(channel, 1)
    for i, b in enumerate(binfactors):
        nsegs = data.size//segmentlength*b

        temp = data[:nsegs*segmentlength//b].reshape(nsegs, segmentlength//b)
        if timeshift > 0:
            t_k1_1 = data[:, :-timeshift].mean(axis=1)
            t_k1_2 = data[:, timeshift:].mean(axis=1)
            t_k2 = (data[:, :-timeshift]*data[:, timeshift:]).mean(axis=1) - t_k1_1*t_k1_2

            q = t_k2/(t_k1_1 + t_k1_2)*2.
            msq[i] = q.mean()
            msq_err[i] = q.std()

        elif timeshift == 0:
            t_k1 = data.mean(axis=1)
            t_k2 = np.power(data, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)

            q = t_k2/t_k1
            msq[i] = q.mean()
            msq_err[i] = q.std()

        else:
            raise ValueError("timeshift is an integer of 0 or positive.")
    return (msq, msq_err)



class MSQTransformer(FFSTransformer):
    """Calculate mean of segmented Q values

    """
    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 8  # default FFS data chunk size for calculating moments

    def __init__(self, segmentlength=DATACYCLE
                     , binfactors=[]
                     , channels=[]
                     , tshift=1):
        super(MSQTransformer, self).__init__(segmentlength, channels)

        if binfactors == []:
            self._binfactors = 2**np.arange(12)
        elif isinstance(binfactors, list):
            self._binfactors = np.array(binfactors).astype(int)
        elif isinstance(binfactors, np.array):
            self._binfactors = binfactors.astype(int)
        else:
            raise TypeError("the type of binfactors is incorrect.")
        if not isinstance(tshift, int):
            raise TypeError("tshift should be an integer value.")
        self._tshift = tshift

    # def fit(self, X, y=None):
    #     return self
    def dark(self, X):
        if not isinstance(X, tmdata):
            raise TypeError("X for dark counts should be a type of tmdata.")
        return dark(X, self.binfactors, self.segmentlength, self.channels[0], self.timeshift)


    def transform(self, X, dark=None):
        """ Calculate MSQ estimators where dark counts can be taken into accout

        dark >> [k1, k2] from dark counts
            (k1 - the first cumulants, k2 - the second cumulants , where each lenghth
            should be the same to the lenghth of binfactors)

        """
        if not isinstance(X, tmdata):
            raise TypeError("X for counts should be a type of tmdata.")
        if dark:
            return msq_est_dark(X, self.binfactors, self.segmentlength, self.channels[0], self.timeshift, dark)
        else:
            return msq_est(X, self.binfactors, self.segmentlength, self.channels[0], self.timeshift)


        # if not isinstance(X, tmdata):
        #     raise TypeError()
        # if self._channels == []:
        #     data = X.data
        #     channels = X.channels
        #     frequency = X.frequency
        #     result = self.calqestimators(X, channels, frequency)
        # else:
        #     data = X.data
        #     channels = X.channels
        #     frequency = X.frequency
        #     result = self.calqestimators(X, self._channels, frequency)
        # return result

    # def calqestimators(self, data, channels, frequency):
    #     result = {}
    #     for i in channels:
    #         if self._tshift == 0:
    #             result[i] = self.sub_qestimators(data[i - 1], frequency)
    #         else:
    #             result[i] = self.sub_qestimators_withtshift(data[i - 1], \
    #                                                 frequency)
    #     return result
    #
    # def sub_qestimators(self, array, frequency):
    #     """Calculate qestimators for tshift == 0 where the shot noise term is
    #     considered.
    #
    #
    #     """
    #     time = np.zeros(self._binfactors.size)
    #     qestimator = np.zeros(self._binfactors.size)
    #     qestimator_stderr = np.zeros(self._binfactors.size)
    #     for i, binf in enumerate(self._binfactors):
    #         time[i] = (self._segmentlength // binf)/float(frequency)
    #
    #         n_segments = array.size // self._segmentlength * binf
    #         temp_data = array[0:n_segments*(self._segmentlength//binf)] \
    #                     .reshape((n_segments, self._segmentlength//binf))
    #         temp_k1 = temp_data.mean(axis=1)
    #         temp_k2 = (temp_data*temp_data).mean(axis=1)
    #         temp_q = temp_k2/temp_k1 - temp_k1 - 1.
    #         qestimator[i] = temp_q.mean()
    #         if n_segments > 1:
    #             qestimator_stderr[i] = temp_q.std()/np.sqrt(n_segments - 1.)
    #         else:
    #             qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1 #TO DO
    #
    #     qestimator_tuple = collections.namedtuple("Qestimator", \
    #                 ["time", "qestimator", "qestimator_stderr", "tsfhit"])
    #     return qestimator_tuple(time, qestimator, qestimator_stderr, self._tshift)
    #
    # def sub_qestimators_withtshift(self, array, frequency):
    #     """Calculate qestimators for tshift >= 1 where no shot noise term exist.
    #
    #     """
    #     time = np.zeros(self._binfactors.size)
    #     qestimator = np.zeros(self._binfactors.size)
    #     qestimator_stderr = np.zeros(self._binfactors.size)
    #     index = 0
    #     for i, binf in enumerate(self._binfactors):
    #         time[i] = (self._segmentlength // binf)/float(frequency)
    #
    #         n_segments = array.size // self._segmentlength * binf
    #         temp_data = (array[0:n_segments*(self._segmentlength//binf)]
    #                     .reshape((n_segments, self._segmentlength//binf)))
    #
    #         mean1 = temp_data[:, :-self._tshift].mean(axis=1).reshape((-1, 1))
    #         mean2 = temp_data[:, self._tshift:].mean(axis=1).reshape((-1, 1))
    #         temp_q = (((temp_data[:, :-self._tshift] - mean1)
    #                     *(temp_data[:, self._tshift:] - mean2)).mean(axis=1)
    #                     .reshape((-1, 1))
    #                     /(mean1 + mean2) * 2 )
    #         qestimator[i] = temp_q.mean()
    #         if n_segments > 1:
    #             qestimator_stderr[i] = temp_q.std()/np.sqrt(n_segments)
    #         else:
    #             qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1
    #
    #     qestimator_tuple = collections.namedtuple("Qestimator", \
    #                 ["time", "qestimator", "qestimator_stderr", "tshift"])
    #     return qestimator_tuple(time, qestimator, qestimator_stderr, self._tshift)


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
