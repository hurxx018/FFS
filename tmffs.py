import numpy as np
from tmdata import tmdata

class tmffs(object):
    """elementary analyses of time-mode photon count data in FFS
    moments, central moments, time-shifted moments, and time-shifted central
    moments

    Q0, Q1
    G0, G1

    """
    CYCLE = 32768   # default FFS data chunk size for saving data from old acquistion card
    DATACYCLE = CYCLE * 4       # default FFS data chunk size for calculating moments

    def __init__(self, ffsdata):
        # check if ffsdata is an instance derived from the superclass tmdata
        if not isinstance(ffsdata, tmdata):
            raise TypeError("ffsdata must be an instance derived from tmdata class")
        self._ffsdata = ffsdata

    def _checkchannel(self, ch):
        """ raise SyntaxError if ch (channel) is not within the measured data channels"""
        if ch not in self.ffsdata.channels:  # check if channel ch within the allowed range
            raise SyntaxError("channel {} NOT in {}".format(ch, self.ffsdata.channels))

    def _checkorder(self, o):
        """ raise SyntaxError if o (order) is not a number"""
        if not isinstance(o, (int, float)):
            raise SyntaxError("order {} has to be a number".format(o))

    def _checkbinfactor(self, lbin):
        """ raise SyntaxError if lbin (length of bin) is not a positive integer"""
        if not (isinstance(lbin, int) and (lbin >= 1)):
            raise SyntaxError("binning factor {} has to be a positive integer".format(lbin))

    def _checktimeshift(self, timeshift):
        """ raise SyntaxError if timeshift is not a non-negative integer"""
        if not (isinstance(timeshift, int) and (timeshift >= 0)):
            raise SyntaxError("timeshift {} has to be a non-negative integer".format(timeshift))

    def _getColumnsRows(self, binneddatasize, nseg, tseg, lseg, lbin):
        """ determine number of columns and rows of binned np.array of length binneddatasize
            Returns: (the number of rows, the number of columns)
            nrowseg : the number of rows as a segment
            ncolseg : the number of columns in a segment
        """
        if isinstance(nseg, int):
            nrowseg = nseg
            ncolseg = binneddatasize//nseg
        elif isinstance(tseg, (int, float)):
            ncolseg = tseg * self._ffsdata.frequency / lbin
            ncolseg = int(ncolseg)
            nrowseg = binneddatasize // ncolseg
        elif isinstance(lseg, int):
            nrowseg = binneddatasize // lseg
            ncolseg = lseg
        else:
            raise SyntaxError("Error in specifying number of segments, length of a segment, time (in s) of a segment: {}, {}, {}".format(nseg,lseg,tseg))
        if nrowseg == 0:
            raise ValueError("Data size is insufficient to reshape the array. Calculated number of rows is zero")
        return (nrowseg, ncolseg)

    def _checkparas(self, ch, o, lbin, timeshift, maxtimeshift, nseg, tseg, lseg):
        """ Check input parameter list:
            ch:     data channel = {1,2,3}
            o:      order (1,2,3,...)
            lbin:   rebin data by factor lbin (positive integer)
            lseg:   segment length of rebinned data (integer)
            tseg:   time period of a segment (in seconds) (float)
            nseg:   number of segments in the data (int)

            """
        self._checkchannel(ch)
        self._checkorder(o)
        self._checkbinfactor(lbin)
        self._checktimeshift(timeshift)

        if maxtimeshift is not None:
            if not isinstance(maxtimeshift, int):
                raise SyntaxError("maxtimeshift {} has to be a non-negative integer".format(maxtimeshift))
            if maxtimeshift < timeshift:
                raise ValueError("maxtimeshift has to be >= timeshift")

        if nseg is not None:
            if not (isinstance(nseg, int) and (nseg >= 1)):
                raise SyntaxError("Number of segment (nseg) {} has to be a positive integer".format(nseg))
        if lseg is not None:
            if not (isinstance(lseg, int) and (lseg >= 1)):
                raise SyntaxError("Length of segment (lseg) {} has to be a positive integer".format(lseg))
        if tseg is not None:
            if type(tseg) in {int, float}:
                if tseg <= 0:
                    raise ValueError("Time period of segment (tseg) {} has to be a positive number".format(tseg))
            else:
                raise SyntaxError("time of segment (tseg) {} has to be a positive integer or a positive float".format(tseg))

    def k(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):
        """ Return the photon counts k^o of channel ch raised by the order o.
            Setting the parameter lbin rebins the photon counts.
                Rebinning = grouping the data in bins of length lbin and summing up the their counts

            Examples:
                k():                return photon counts k of channel 1
                k(o=2):             return squared photon counts k^2 of channel 1:
                k(ch=1, lbin=10):   return photon counts of channel 1 rebinned by a factor 10

            In addition, the (rebinned) data
                """
        self._checkparas(ch=ch, o=o, lbin=lbin, timeshift=timeshift, maxtimeshift=maxtimeshift, nseg=nseg, tseg=tseg, lseg=lseg)
        binneddatasize = self._ffsdata.data[ch - 1].size//lbin
        # get photon counts k1 of channel ch
        if lbin > 1:
            # data are first rebiined by a factor lbin
            k1 = self._ffsdata.data[ch-1][0:binneddatasize*lbin].reshape(binneddatasize, lbin).sum(axis=1)
        else:
            k1 = self._ffsdata.data[ch-1]
        if timeshift > 0:
            k1 = k1[timeshift:]
        if maxtimeshift is not None:
                nk1 = k1.size
                k1 = k1[:nk1-maxtimeshift+timeshift]
        binneddatasize = k1.size
        if not flatten:
            # if lseg, tseg or nseg present, then reshape the data into a matrix k1[nseg, lseg)
            if not (lseg is None and tseg is None and nseg is None):
                nrow, ncol = self._getColumnsRows(binneddatasize, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
                k1 = k1[0:nrow*ncol].reshape(nrow, ncol)
        # finally return k1^o
        if o == 1:
            return k1
        else:
            return np.power(k1, o)

    def time(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):
        self._checkparas(ch=ch, o=o, lbin=lbin, timeshift=timeshift, maxtimeshift=maxtimeshift, nseg=nseg, tseg=tseg, lseg=lseg)
        Ts = lbin/self.ffsdata.frequency               # binned sampling time
        binneddatasize = self._ffsdata.data[ch-1].size//lbin
        time = np.arange(binneddatasize)*Ts
        if timeshift > 0:
            time = time[timeshift:]
        if maxtimeshift is not None:
                time = time[:binneddatasize-maxtimeshift+timeshift]
        binneddatasize = time.size
        if not flatten:
            if not (lseg is None and tseg is None and nseg is None):
                nrow, ncol = self._getColumnsRows(binneddatasize, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
                time = time[0:nrow*ncol].reshape(nrow, ncol)
        if o == 1:
            return time
        else:
            return np.power(time, o)

    def mtime(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE):
        return self.time(ch, o, lbin, timeshift, maxtimeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def mk(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE):
        #nsegCol, lsegRow = self._checkparas(ch, o, lbin, nseg, tseg, lseg)
        return self.k(ch, o, lbin, timeshift, maxtimeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)
        """ where time-shifted """

    def dk(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):
        """calculate the fluctuations    ;;; ordered
        """
        k1  = self.k(ch, o=1, lbin=lbin, timeshift=timeshift, maxtimeshift=maxtimeshift, nseg=nseg, tseg=tseg, lseg=lseg, flatten=False)
        mk1 = self.mk(ch, o=1, lbin=lbin, timeshift=timeshift, maxtimeshift=maxtimeshift, nseg=nseg, tseg=tseg, lseg=lseg)
        dk1 = k1 - mk1[:, np.newaxis]
        if o != 1:
            dk1 = np.power(dk1, o)
        if flatten:
            return dk1.flatten()
        else:
            return dk1

    def mdk(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE):
        return self.dk(ch, o, lbin, timeshift, maxtimeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def kk(self, ch=(1, 1), o=(1, 1), lbin=1, timeshift=(0, 0), nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):

        maxts = max(timeshift)
        k0 = self.k(ch=ch[0], o=o[0], lbin=lbin, timeshift=timeshift[0], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        k1 = self.k(ch=ch[1], o=o[1], lbin=lbin, timeshift=timeshift[1], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        # check here
        k0k1 = k0*k1
        #k0k1 = k0[:k0.size-timeshift]*k1[timeshift:]
        if flatten:
            return k0k1
        else:
            nrow, ncol = self._getColumnsRows(k0k1.size, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
            return k0k1[0:nrow*ncol].reshape(nrow, ncol)

    def mkk(self, ch=(1, 1), o=(1, 1), lbin=1, timeshift=(0, 0), nseg=None, tseg=None, lseg=DATACYCLE):
        return self.kk(ch, o, lbin, timeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def dkdk(self, ch=(1, 1), o=(1, 1), lbin=1, timeshift=(0, 0), nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):

        maxts = max(timeshift)
        dk0 = self.dk(ch=ch[0], o=o[0], lbin=lbin, timeshift=timeshift[0], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        dk1 = self.dk(ch=ch[1], o=o[1], lbin=lbin, timeshift=timeshift[1], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        dk0dk1 = dk0*dk1
        # dk0dk1 = dk0[:dk0.size-timeshift]*dk1[timeshift:]  # this one should be corrected: The calculation of mean value should be affected by the time shift factor.
        if flatten:
            return dk0dk1
        else:
            nrow, ncol = self._getColumnsRows(dk0dk1.size, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
            return dk0dk1[0:nrow*ncol].reshape(nrow, ncol)

    def mdkdk(self, ch=(1, 1), o=(1, 1), lbin=1, timeshift=(0, 0), nseg=None, tseg=None, lseg=DATACYCLE):
        return self.dkdk(ch, o, lbin, timeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def kkk(self, ch=(1, 1, 1), o=(1, 1, 1), lbin=1, timeshift=(0, 0, 0), nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):

        maxts = max(timeshift)
        k0 = self.k(ch=ch[0], o=o[0], lbin=lbin, timeshift=timeshift[0], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        k1 = self.k(ch=ch[1], o=o[1], lbin=lbin, timeshift=timeshift[1], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        k2 = self.k(ch=ch[2], o=o[2], lbin=lbin, timeshift=timeshift[2], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        #nk0 = k0.size
        #nk1 = k1.size
        #nk2 = k2.size
        #check that nk0 == nk1 == nk2
        #maxTS = max(timeshift)
        kkk = k0*k1*k2
        if flatten:
            return kkk
        else:
            nrow, ncol = self._getColumnsRows(kkk.size, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
            return kkk[0:nrow*ncol].reshape(nrow, ncol)

    def mkkk(self, ch=(1, 1, 1), o=(1, 1, 1), lbin=1, timeshift=(0, 0, 0), nseg=None, tseg=None, lseg=DATACYCLE):
        return self.kkk(ch, o, lbin, timeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def dkdkdk(self, ch=(1, 1, 1), o=(1, 1, 1), lbin=1, timeshift=(0, 0, 0), nseg=None, tseg=None, lseg=DATACYCLE, flatten=True):

        maxts = max(timeshift)
        dk0 = self.dk(ch=ch[0], o=o[0], lbin=lbin, timeshift=timeshift[0], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        dk1 = self.dk(ch=ch[1], o=o[1], lbin=lbin, timeshift=timeshift[1], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        dk2 = self.dk(ch=ch[2], o=o[2], lbin=lbin, timeshift=timeshift[2], maxtimeshift=maxts, nseg=nseg, tseg=tseg, lseg=lseg)
        #nk0 = dk0.size
        #nk1 = dk1.size
        #nk2 = dk2.size
        #check that nk0 == nk1 == nk2
        #maxTS = max(timeshift)
        dkkk = dk0*dk1*dk2
        if flatten:
            return dkkk
        else:
            nrow, ncol = self._getColumnsRows(dkkk.size, nseg=nseg, tseg=tseg, lseg=lseg, lbin=lbin)
            return dkkk[0:nrow*ncol].reshape(nrow, ncol)

    def mdkdkdk(self, ch=(1, 1, 1), o=(1, 1, 1), lbin=1, timeshift=(0, 0, 0), nseg=None, tseg=None, lseg=DATACYCLE):
        return self.dkdkdk(ch, o, lbin, timeshift, nseg, tseg, lseg, flatten=False).mean(axis=1)

    def Q0(self, ch=1, lbin=1, nseg=None, tseg=None, lseg=DATACYCLE):
        k1 = self.mk(ch=ch, o=1, lbin=lbin, timeshift=0, nseg=nseg, tseg=tseg, lseg=lseg)
        dk2 = self.mdk(ch=ch, o=2, lbin=lbin, timeshift=0, nseg=nseg, tseg=tseg, lseg=lseg)
        return dk2/k1 - 1

    def Q1(self, ch=1, lbin=1, nseg=None, tseg=None, lseg=DATACYCLE):
        k1a = self.mk(ch=ch, o=1, lbin=lbin, timeshift=0, maxtimeshift=1, nseg=nseg, tseg=tseg, lseg=lseg)
        k1b = self.mk(ch=ch, o=1, lbin=lbin, timeshift=1, maxtimeshift=1, nseg=nseg, tseg=tseg, lseg=lseg)
        dk2 = self.mdkdk(ch=(ch, ch), o=(1, 1), lbin=lbin, timeshift=(0,1), nseg=nseg, tseg=tseg, lseg=lseg)
        return dk2/np.sqrt(k1a*k1b)

    def g0(self, ch=1, lbin=1, nseg=None, tseg=None, lseg=DATACYCLE):
        k1 = self.mk(ch=ch, o=1, lbin=lbin, timeshift=0, nseg=nseg, tseg=tseg, lseg=lseg)
        Q0 = self.Q0(ch=ch, lbin=lbin, nseg=nseg, tseg=tseg, lseg=lseg)
        return Q0/k1

    def g1(self, ch=1, lbin=1, nseg=None, tseg=None, lseg=DATACYCLE):
        k1a = self.mk(ch=ch, o=1, lbin=lbin, timeshift=0, maxtimeshift=1, nseg=nseg, tseg=tseg, lseg=lseg)
        k1b = self.mk(ch=ch, o=1, lbin=lbin, timeshift=1, maxtimeshift=1, nseg=nseg, tseg=tseg, lseg=lseg)
        dk2 = self.mdkdk(ch=(ch, ch), o=(1, 1), lbin=lbin, timeshift=(0,1), nseg=nseg, tseg=tseg, lseg=lseg)
        return dk2/(k1a*k1b)

    def mF(self, ch=1, o=1, lbin=1, timeshift=0, maxtimeshift=None, nseg=None, tseg=None, lseg=DATACYCLE):
        return self.mk(ch, o, lbin, timeshift, maxtimeshift, nseg, tseg, lseg)*self._ffsdata.frequency/lbin

    @property
    def ffsdata(self):
        return self._ffsdata


def main():
    import matplotlib.pyplot as plt
    from readFFSfrombinfiles import readFFSfrombinfiles

    filename1 = "A488_cal.1.001.bin"
    filename2 = "A488_cal.2.001.bin"
    filename3 = "A488_cal.3.001.bin"
    data = readFFSfrombinfiles([filename1, filename2, filename3], [1, 2, 3], frequency=100000)
    #data._data = [np.array([1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2]), np.array([2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1]),np.arange(28)]
    # data._data = [np.arange(20)]
    d = tmffs(data)

    dkk1 = d.dkdk(timeshift=(0,1),nseg=2)
    kkk1 = d.kkk(timeshift=(0,1,2))
    kkk2 = d.kkk(timeshift=(0,0,0))
    dkkk1 = d.dkdkdk(timeshift=(0,1,2),nseg=2)
    time1 = d.time(timeshift=1)
    time1sqr = d.time(o=2)
    # rebin photoncounts of ch 1 by 4 [K, K, ...] with K = k+k+k+k
    k1 = d.k(lbin=4)
    print(k1)
    # square of photoncounts [k^2, k^2, k^2, ...] from ch 1
    k2 = d.k(o=2)
    # photoncounts from ch 2
    k3 = d.k(ch=2)
    # photoncounts k divided into 2 segments
    k4 = d.k(nseg=4, flatten=False)
    print(k4)
    # note that segmenting of the photoncounts is used primarily for internal calculations
    # thus, the method flattens the array before returning it
    # if you want to see the segmenting, you have to turn flatten off
    k5 = d.k(nseg=2, flatten=False)
    # the above returns a matrix[nrow, ncol] with nrow containing each segment for a total of ncol segments
    # let's move on to the mean of the photon counts
    # mean photon counts per segment (by default a segment = 32768*4 data points) of ch 1
    mk1 = d.mk()
    # mean photoncounts <k> of channel 2 over segments of 2 seconds
    mk2 = d.mk(ch=2, tseg=2)
    # mean rebinned photoncounts over segments of 1 second. Binfactor is 4 and channel is 1
    mk3 = d.mk(lbin=4, tseg=1)
    # mean squared photoncounts of segments of length 10000. Channel = 2
    mk4 = d.mk(ch=2, o=2, lseg=1000000)
    # now calculate fluctuations
    # calculate dk = k - <k>_seg. data are segmented, the average of the segment is subtracted from individual photoncounts k. data are from channel 1, segment length is 2 second.
    dk1 = d.dk(tseg=2)
    # calculate dk^2
    dk2 = d.dk(o=2, tseg=2)
    # bin data of channel 2 by a factor of 10, then calculate dk^2 using
    dk3 = d.dk(o=2, lbin=10, nseg=2)
    # calculate the variance of the photon counts from channel 3 reboinned by a factor of 5 using a segment time of 0.5 seconds
    mdk1 = d.mdk(o=3, lbin=5, tseg=0.5)
    # make time vectors
    # time vector of data
    # time vector for channel 1 in units of seconds ([0, Ts, ...] with Ts = sampling time)
    t1 = d.time()
    # time vector in units of second for rebinned data ([0, 10Ts, 20Ts, ...] with 10Ts = Tb as binning time)
    t2 = d.time(lbin = 10)
    # averaged time over two segments  ([<time>_seg, <time>_seg])
    mt1 = d.mtime(nseg=2)
    # plot average photon counts per sample versus segment time
    plt.plot(d.mtime(tseg=1), d.mk(tseg=1))
    plt.show()
    # determine  k1[i]*k1[i], k1[i+1]*k1[i+1] of channel 2
    # which is the same as d.k(ch=1, o=2)
    kk1 = d.kk(ch=(1, 1))
    # or as moment <k1 k1> = <k1^2> per segment
    mkk1 = d.mkk(ch=(1, 1), tseg=1)
    # determine the bivariate moment (1,1) between the 1st and 2nd channel
    mkk2 = d.mkk(ch=(1, 2), tseg=1)
    # determine the bivariate moment (1,2) between the 1st and 3rd channel [<k1 k3^2>_seg, <k1 k3^2>_seg, ...]
    mkk3 = d.mkk(ch=(1,3), o=(1, 2), tseg=1)
    # calculate timeshifted moment <k1[i] k1[i+1]> per segment
    mkk3 = d.mkk(timeshift=(0,1), tseg=1)
    # calculate timeshifted moment <k1[i] k2[i+1]^2> per segment
    mkk4 = d.mkk(ch=(1, 2), o=(1, 2), timeshift=(0,1), tseg=1)
    # calculate timeshifted moment <K1[i] K2[i+1]^2> per segment for data rebinned by a factor of 5 (K1[i] = k1+k1+k1+k1+k1)
    mkk5 = d.mkk(ch=(1, 2), o=(1, 2), timeshift=(0,1), tseg=1, lbin=5)
    #similarily, we can calculate 2-order fluctuations
    # dk1[i]^2*dk2[i]^2
    dkdk1 = d.dkdk(ch=(1, 2),o=(2, 2), tseg=1)
    # or timeshifted  dk1[i]dk1[i+1]
    dkdk2 = d.dkdk(ch=(1, 1),o=(1, 1), timeshift=(0,1), tseg=1)
    # second central moments
    # <dk1[i] dk2[i+1]>_seg
    mdkdk1 = d.mdkdk(ch=(1,2),timeshift=(0,1), tseg=4)
    #third order terms and third order moments
    # k1[i] k2[i] k3[i+1]
    kkk1 = d.kkk(ch=(1, 2, 3), timeshift=(0, 0, 1))
    # <k1[i] k2[i] k3[i+1]>_seg
    mkkk1 = d.mkkk(ch=(1, 2, 3), timeshift=(0, 0, 1), tseg=1)
    #third order (central) terms and third order central moments
    # dk1[i] dk2[i] dk3[i+1]
    dkkk1 = d.dkdkdk(ch=(1, 2, 3), timeshift=(0, 0, 1))
    # <dk1[i] dk2[i] dk3[i+1]>_seg
    mdkkk1 = d.mdkdkdk(ch=(1, 2, 3), timeshift=(0, 0, 1), tseg=1)
    # Q0 = Mandel's Q parameter
    # Q of channel 2 of rebinned data
    Q0values = d.Q0(ch=2, lbin=5, tseg=1)
    # timeshifted Q of channel 2 of rebinned data
    Q1values = d.Q1(ch=2, lbin=5, tseg=1)
    # g0 of channel 2
    g0values = d.g0(ch=2, tseg=5)
    # g1 of channel 2
    g1values = d.g1(ch=2, tseg=5)
    # test:
    # return the average intenity per segment (counts per second)
    F1 = d.mF(tseg=1)
    print(F1)


if __name__ == "__main__":
    main()
