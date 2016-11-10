import numpy as np
from zscan_tools/zscanPSFGLsp import zscanPSFGLsp



class zscanFitter(object):
    """docstring for zscanFitter."""

    #TO DO: explain each geo in detail
    geodict = {"BASE":-1, "SLAB(L)":0, "UP":1, "DOWN":2, "DELTA":3
                , "SPHERE(R)":4, "CYLINDER(R)":5
                , "notSLAB":10, "notDELTA":13, "notSPHERE":14, "notCYLINDER":15
                , "notSPHERE:up":24, "notCYLINDERup":25}
    # mGL : zR, y, w0
    # 3DG : zR, w0
    # GL  : zR, w0
    psfdict = {"mGL":[], "3DG":[], "GL":[]}



    def __init__(self):
        super(zscanFitter, self).__init__()
        self._psfmodel = ""     # psfmodel
        self._psf = {}          # for psf paras on each channel
        self._psffixed = {}     # fixed condition for psf
        self._geoinfo = {}      # for sample geometric models
        self._geoinfofixed = {}

    def setkz(self, z, kz=[], channels=[]):
        # TODO
        # self.plot()
        pass

    def fit(self, z, kz=[], channels=[]):
        # TODO build a fit procedure
        pass


    def setPSF(self, channel=1, psfparas=[], fixed=[], psfmodel="mGL"):
        """set the PSF model with psf paras and fixed conditions

        psfparas = [zR, y, w0] for mGL
        psfparas = [zR, w0] for GL
        psfparas = [zR, w0] for 3DG

        """
        if self._psfmodel == "":
            # initialize the psf model
            self._psfmodel = psfmodel
        elif self._psfmodel != psfmodel:
            raise ValueError("psfmodel is not consistent.")

        if self._psfmodel == "mGL":
            if psfparas != []:
                assert (len(psfparas) ==  3), "The number of paras should be 3."
        elif self._psfmodel == "GL":
            if psfparas != []:
                assert (len(psfparas) ==  2), "The number of paras should be 2."
        elif self._psfmodel == "3DG":
            if psfparas != []:
                assert (len(psfparas) ==  2), "The number of paras should be 2."

        if fixed == []:
            fixed = [0]*len(psfparas)
            #fixed = [False]*len(psfparas)
        else:
            assert (len(paras) ==  len(fixed)), \
                "The number of paras does not match to the number of fixed."

        self._psf[channel] = psfparas
        self._psffixed[channel] = fixed
        return

    def getPSF(self, channel=1):
        return self._psfmodel, self._psf[channel], self._psffixed[channel]

    def addLayer(self, geomodel, channel=1):
        """geomodel is either integer or string in geodict."""
        tempmodel = geomodel
        if isinstance(tempmodel, int):
            for key, value in self.geodict.items():
                if value == tempmodel:
                    tempmodel = key
                    break
        self.setLayer(tempmodel, channel=channel)
        return

    def setLayer(self, geomodel, channel=1, layer_index=None,
                            F=0., b=0., LR=0., fF=0, fb=0, fLR=0):
        """set geo model on each layer with geo-paras and fixed conditions

        F >> counts per bin
        b >> ??
        LR >> length or radius

        fF >> 0 (free) or 1 (fixed) for F
        b >> 0 (free) or 1 (fixed) for b
        fLR >> 0 (free) or 1 (fixed) for LR
        """
        if layer_index == None:
            self._geoinfo[channel].append([geomodel, F, b, LR])
            self._geoinfofixed[channel].append([geomodel, fF, fb, fLR])
        elif (isinstance(layer_index, int) and
                layer_index < len(self._geoinfo[channel])):
            if self._geoinfo[channel][layer_index][0] == geomodel:
                self._geoinfo[channel][layer_index][1:4] = [F, b, LR]
                self._geoinfofixed[channel][layer_index][1:4] = [fF, fb, fLR]
            else:
                raise ValueError("model does not match.")
        else:
            raise ValueError("layer_index is out of the range")
        return


    def removeLayer(self, layer_index=None):
        # remove a layer from all channels
        for key in self._geoinfo:
            if self._geoinfo[key] != []:
                try:
                    if isinstance(layer_index, int):
                        self._geoinfo[key].pop(layer_index)
                        self._geoinfofixed[key].pop(layer_index)
                    elif layer_index == None:
                        self._geoinfo[key].pop()
                        self._geoinfofixed[key].pop()
                    else:
                        raise ValueError("layer_index is out of the allowed range")
                except:
                    raise ValueError("layer_index is out of the allowed range")
            else:
                print("geoinfo is empty.")
        return


    def fitkzMultiLayerFCT(cls, z, paras):
        pass


    @property
    def psfmodel(self):
        return self._psfmodel
    @property
    def psf(self):
        return self._psf
    @property
    def psffixed(self):
        return self._psffixed
    @property
    def geoinfo(self):
        return self._geoinfo
    @property
    def geoinfofixed(self):
        return self._geoinfofixed

    @classmethod
    def printgeodict(cls):
        temp = sorted(cls.geodict.items(), key= (lambda x:x[1]))
        for key, value in temp:
            print(key, value)
