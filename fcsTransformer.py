import itertools

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
    def __init__(self, segmentlength=32768L*4, channels=[]):
        self._segmentlength = segmentlength
        self._channels = channels

    def fit(self, X, y=None):
        return self

    def transform(self, X):

        if self._channels == []:
            X.get_data()
            X.get_channels()

        else:
        #
            pass


        return result

    def cal_correlations(data, channels, segmentlength):
        n_segment = data[0].size//segmentlength
        result = {}
        for x in itertools.combinations_with_replacement(channels, 2):
            arr1 = data[x[0]][0:segmentlength*n_segment].reshape((n_segment, segmentlength))
            arr2 = data[x[1]][0:segmentlength*n_segment].reshape((n_segment, segmentlength))
            result[x] = sub_correlations(arr1, arr2, segmentlength)
        return result

    def sub_correlations(arr1, arr2, segmentlength):


        x1 = arr1.mean(axis=1)
        np.repeat(x1)



    def rebin_factor( a, newshape ):
        '''Rebin an array to a new shape.
        newshape must be a factor of a.shape.
        '''
        assert len(a.shape) == len(newshape)
        assert not sometrue(mod( a.shape, newshape ))

        slices = [ slice(None,None, old/new) for old,new in zip(a.shape,newshape) ]
        return a[slices]
