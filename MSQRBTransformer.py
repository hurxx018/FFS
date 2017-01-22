from collections import namedtuple
from itertools import combinations_with_replacement as cwr
from itertools import product

import numpy as np

from tmdata import tmdata, tmdataSR, tmdataMR
from FFSTransformer import FFSTransformer

class MSQRBTransformer(FFSTransformer):
    """Calculate mean of segmented Q values with Rebinning

    """
    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 8  # default FFS data chunk size for calculating moments

    def __init__(self, segmentlength=DATACYCLE
                     , binfactors=[]
                     , channels=[]
                     , rbfs=[]  # rebin-factors
                     , tshift=[]):
        super(MSQRBTransformer, self).__init__(segmentlength, channels)

        if not binfactors:
            self._binfactors = 2**np.arange(12)
        elif isinstance(binfactors, list):
            self._binfactors = np.array(binfactors).astype(int)
        elif isinstance(binfactors, np.array):
            self._binfactors = binfactors.astype(int)
        else:
            raise TypeError("the type of binfactors is incorrect.")

        if not rbfs:
            self._rbfs = np.array([1,2,4,8,16])
        elif isinstance(rbfs, list):
            self._rbfs = np.array(rbfs, dtype=int)
        else:
            raise TypeError("the type of rbfs is incorrect.")

        if not tshitf:
            self._tshift = np.array([1])
        elif isinstance(tshift, list):
            self._tshift = np.array(tshift, dtype=int)
        else:
            raise TypeError("The type of tshift is incorrect.")

    # def fit(self, X, y=None):
    #     return self

    def transform(self, X):
        """Calculate MSQ estimators from input X with Rebinnings.

        """
        if isinstance(X, (tmdata, tmdataSR)):
            temp_X = self._gettmdataMR(X)
        elif isinstance(X, tmdataMR):
            temp_X = X
        else:
            raise TypeError("The type of input X should be one of tmdata"
                + ", tmdataSR, and tmdataMR")

        if self._channels == []:
            pairs_of_cha = cwr(temp_X.channels, 2)  # all pairs of channel numbers
        else:
            pairs_of_cha = cwr(self._channels, 2)
        return self.calqestimators(temp_X, pairs_of_cha, temp_X.frequency)
    #
    def calqestimators(self, data, pairs, frequency):
        tvec, lvec = None, None
        tvec = np.zeros(self._binfactors.size)
        for i, binf in enumerate(self._binfactors):
            tvec[i] = (self._segmentlength // binf)/float(frequency)

        result = {}
        for rbf, tshift in product(self._rbfs, self._tshift):
            for pair in pairs:
                k0 = data.rebin(cha=pair[0], lbin=rbf)
                k1 = data.rebin(pair[1], rbf)
                key = (rbf, tshift, pair[0], pair[1])
                if tshift == 0:
                    result[key] = self._sub_qestimators(k0, k1, rbf, pair)
                else:
                    result[key] = self._sub_qestimators_withtshift(k0, k1, \
                                                        rbf, tshift)
        return result

    def _sub_qestimators(self, arr0, arr1, rbf, pair):
        """Calculate qestimators for tshift == 0 where the shot noise term is
        considered if arr0 == arr1.
    #
    #
    #     """
        qestimator = np.zeros(self._binfactors.size)
        qestimator_stderr = np.zeros(self._binfactors.size)
        for i, binf in enumerate(self._binfactors):
            new_sgl = self._segmentlength // rbf // binf
            n_segments = arr0.size // new_sgl

            td0 = (arr0[0:n_segments*new_sgl]
                        .reshape((n_segments, new_sgl)))
            td1 = (arr1[0:n_segments*new_sgl]
                        .reshape((n_segments, new_sgl)))

            mean0 = td0.mean(axis=1).reshape((-1, 1))
            mean1 = td1.mean(axis=1).reshape((-1, 1))

            if pair[0] ==  pair[1]:
                temp_q = (((td0 - mean0)*(td1 - mean1)).mean(axis=1)
                            .reshape((-1, 1)) / mean0) - 1.
            else:
                temp_q = (((td0 - mean0)*(td1 - mean1)).mean(axis=1)
                            .reshape((-1, 1)) / mean0)

            qestimator[i] = temp_q.mean()
            if n_segments > 1:
                qestimator_stderr[i] = temp_q.std()/np.sqrt(n_segments - 1.)
            else:
                qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1 #TO DO
    #
        qestimator_tuple = collections.namedtuple("Qestimator", \
                    ["time", "qestimator", "qestimator_stderr"])
        return qestimator_tuple(time, qestimator, qestimator_stderr)
    #
    def _sub_qestimators_withtshift(self, arr0, arr1, rbf, tshift):
        """Calculate qestimators for tshift >= 1 where no shot noise term exist.

        """
        # time = np.zeros(self._binfactors.size)
        qestimator = np.zeros(self._binfactors.size)
        qestimator_stderr = np.zeros(self._binfactors.size)
        # index = 0
        for i, binf in enumerate(self._binfactors):
            # time[i] = (self._segmentlength // binf)/float(frequency)
            new_sgl = self._segmentlength // rbf // binf
            n_segments = arr0.size // new_sgl

            td0 = (arr0[0:n_segments*new_sgl]
                        .reshape((n_segments, new_sgl)))
            td1 = (arr1[0:n_segments*new_sgl]
                        .reshape((n_segments, new_sgl)))

            mean0 = td0[:, :-tshift].mean(axis=1).reshape((-1, 1))
            mean1 = td1[:, tshift:].mean(axis=1).reshape((-1, 1))

            temp_q = (((td0[:, :-tshift] - mean0)
                        *(td1[:, tshift:] - mean1)).mean(axis=1)
                        .reshape((-1, 1)) / mean0 )
            qestimator[i] = temp_q.mean()
            if n_segments > 1:
                qestimator_stderr[i] = temp_q.std()/np.sqrt(n_segments - 1)
            else:
                qestimator_stderr[i] = np.abs(qestimator[i]) * 0.1

        qestimator_tuple = namedtuple("Qestimator", \
                    ["qestimator", "qestimator_stderr"])
        return qestimator_tuple(qestimator, qestimator_stderr)


    @property
    def binfactors(self):
        return self._binfactors

    @binfactors.setter
    def binfactors(self, values):
        if isinstance(values, np.ndarray):
            self._binfactors = values
        else:
            self._binfactors = np.ndarray(values, dtype=int)

    @property
    def tshift(self):
        return self._tshift

    @tshift.setter
    def tshift(self, value):
        self._tshift = value
