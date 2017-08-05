import collections
import numpy as np

from tmdata import tmdata
from FFSTransformer import FFSTransformer

def dark(X_dark, binfactors, segmentlength=32768*4, channel=1, tshift=1):
    """ Calculate k1 and k2 from tmdata in terms of binfactors

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    tshift >>  a factor of timeshift; an integer

    Output >> (first order, second order) factorial cumulants

    """
    k1 = np.zeros( len(binfactors) )
    k2 = np.zeros( len(binfactors) )
    data = X_dark.ch(channel, 1)
    for i, b in enumerate(binfactors):
        nsegs = data.size//segmentlength*b
        temp = data[:nsegs*segmentlength//b].reshape(nsegs, segmentlength//b)
        if tshift > 0:
            t_k1_1 = temp[:, :-tshift].mean(axis=1)
            t_k1_2 = temp[:, tshift:].mean(axis=1)
            t_k2 = (temp[:, :-tshift]*temp[:, tshift:]).mean(axis=1) - t_k1_1*t_k1_2
            k1[i] = (t_k1_1 + t_k1_2).mean(axis=0)/2.
            k2[i] = t_k2.mean(axis=0)
        elif tshift == 0:
            t_k1 = data.mean(axis=1)
            t_k2 = np.power(data, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)
            k1[i] = t_k1.mean(axis=0)
            k2[i] = t_k2.mean(axis=0)
        else:
            raise ValueError("tshift is an integer of 0 or positive.")
    return (k1, k2)


def msq_est_dark(X, binfactors, segmentlength=32768*4, channel=1, tshift=1, dark=None):
    """ Calculate MSQ estimators with dark count corrections

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    tshift >>  a factor of timeshift; an integer

    dark >> [k1, k2] from a function dark


    Output >> (tseg, msq, msq_err)  ; No consideration of frequency

    """
    if dark:
        assert len(binfactors) == dark[0].size, "inconsistency in the size of dark k1 and k2"
    else:
        raise ValueError("dark is not available.")

    tseg = np.zeros( len(binfactors) )
    msq = np.zeros( len(binfactors) )
    msq_err = np.zeros( len(binfactors) )

    data = X.ch(channel, 1)
    for i, binf in enumerate(binfactors):
        tseg[i] = segmentlength // binf
        nsegs = data.size//segmentlength*binf

        dark1 = np.ones(nsegs) * dark[0][i]
        dark2 = np.ones(nsegs) * dark[1][i]

        temp = data[:nsegs*segmentlength//binf].reshape(nsegs, segmentlength//binf)
        if tshift > 0:
            t_k1_1 = temp[:, :-tshift].mean(axis=1)
            t_k1_2 = temp[:, tshift:].mean(axis=1)
            t_k2 = (temp[:, :-tshift]*temp[:, tshift:]).mean(axis=1) - t_k1_1*t_k1_2

            m1 = (t_k1_1 + t_k1_2)/2. - dark1   # dark count correction
            m2 = t_k2 - dark2                   # dark count correction
            q = m2/m1
            msq[i] = q.mean()
            msq_err[i] = q.std()/np.sqrt(nsegs-1)

        elif tshift == 0:
            t_k1 = temp.mean(axis=1)
            t_k2 = np.power(temp, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)

            m1 = t_k1 - dark1
            m2 = t_k2 - dark2
            q = m2/m1
            msq[i] = q.mean()
            msq_err[i] = q.std()/np.sqrt(nsegs-1)

        else:
            raise ValueError("tshift is an integer of 0 or positive.")
    return (tseg, msq, msq_err)


def msq_est(X, binfactors, segmentlength=32768*4, channel=1, tshift=1):
    """Calculate MSQ estimators without dark count corrections

    X_dark >> tmdata for dark counts
    binfactors >> list
    segmentlength >> default 32768*4
    channel >>
    tshift >>  a factor of timeshift; an integer

    Output >> (tseg, msq, msq_err)  ; No consideration of frequency

    """

    tseg = np.zeros( len(binfactors) )
    msq = np.zeros( len(binfactors) )
    msq_err = np.zeros( len(binfactors) )

    data = X.ch(channel, 1)
    for i, binf in enumerate(binfactors):
        tseg[i] = segmentlength // binf
        nsegs = data.size//segmentlength*binf
        temp = data[ : nsegs*segmentlength//binf]
        temp = temp.reshape(nsegs, segmentlength//binf)
        if tshift > 0:
            t_k1_1 = temp[:, :-tshift].mean(axis=1)
            t_k1_2 = temp[:, tshift:].mean(axis=1)
            t_k2 = (temp[:, :-tshift]*temp[:, tshift:]).mean(axis=1) - t_k1_1*t_k1_2

            q = t_k2/(t_k1_1 + t_k1_2)*2.
            msq[i] = q.mean()
            msq_err[i] = q.std()/np.sqrt(nsegs-1)

        elif tshift == 0:
            t_k1 = data.mean(axis=1)
            t_k2 = np.power(data, 2).mean(axis=1) - t_k1 - np.power(t_k1, 2)

            q = t_k2/t_k1
            msq[i] = q.mean()
            msq_err[i] = q.std()/np.sqrt(nsegs-1)

        else:
            raise ValueError("tshift is an integer of 0 or positive.")
    return (tseg, msq, msq_err)



class MSQTransformer(FFSTransformer):
    """Calculate mean of segmented Q estimators from photon count data

    Currently only for single channel

    segmentlength >> an integer
    binfactors    >> list/np.array
    channels      >> list e.g. [1]
    tshift        >> a factor of time-shift (default : 1) which can be zero or
                    a positive integer

    """
    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 8  # default FFS data chunk size for calculating moments

    def __init__(self, segmentlength=DATACYCLE
                     , binfactors=[]
                     , channels=[]
                     , tshift=1):
        super(MSQTransformer, self).__init__(segmentlength, channels)

        if binfactors == []:
            self._binfactors = 2**np.arange(10, dtype=np.int64)
        elif isinstance(binfactors, list):
            self._binfactors = np.array(binfactors, dtype=np.int64)
        elif isinstance(binfactors, np.array):
            self._binfactors = binfactors.astype(np.int64)
        else:
            raise TypeError("the type of binfactors is incorrect.")
        if not isinstance(tshift, int):
            raise TypeError("tshift should be an integer value.")
        self._tshift = tshift

    # def fit(self, X, y=None):
    #     return self
    def dark(self, X):
        """ Calculate the first and second order factorial cumulants from
        dark-count data

        X >> tmdata including dark counts

        Outputs >> (first order, second order)

        """
        if not isinstance(X, tmdata):
            raise TypeError("X for dark counts should be a type of tmdata.")
        if not self._channels:
            tch = self._channels[0]
        else:
            tch = X.channels[0]

        return dark(X, self._binfactors, self._segmentlength, tch, self._tshift)


    def transform(self, X, dark=None):
        """ Calculate MSQ estimators where dark counts can be taken into accout

        dark >> [k1, k2] from dark counts
            (k1 - the first cumulants, k2 - the second cumulants , where each lenghth
            should be the same to the lenghth of binfactors)

        Output >> a namedtuple (tseg, msq_est, msq_err, tshift)
        """
        if not isinstance(X, tmdata):
            raise TypeError("X for counts should be a type of tmdata.")
        if not self._channels:
            tch = self._channels[0]
        else:
            tch = X.channels[0]

        if dark:
            temp = msq_est_dark(X, self._binfactors, self._segmentlength, tch, self._tshift, dark)
        else:
            temp = msq_est(X, self._binfactors, self._segmentlength, tch, self._tshift)
        msq_est_tuple = collections.namedtuple("MSQ", \
                            ["t", "msq", "err", "tsfhit"])
        return msq_est_tuple(temp[0]/X.frequency, temp[1], temp[2], self._tshift)

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

    from rFFS import rFFS
    import matplotlib.pyplot as plt

    filename2 = "A488_cal.2.001.bin"
    ffs = rFFS( [1], 100000, "iss0" )
    tdata = ffs.load([filename2])
    msq1 = MSQTransformer(segmentlength=32768*8*4, channels=[1], tshift=1)

    msq1_est = msq1.transform(tdata)
    plt.plot(msq1_est[0], msq1_est[1], '*')
    plt.xscale('log')
    plt.show()

    dark = msq1.dark(tdata)
    msq1_est = msq1.transform(tdata, dark)
    plt.plot(msq1_est[0], msq1_est[1], '*')
    plt.xscale('log')
    plt.show()

if __name__ == "__main__":
    main()
