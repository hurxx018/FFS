import numpy as np
from zscan_tools.zscanPSFGLsp import zscanPSFGLsp
from zscan_tools.zscanMultiLayer import zscanMultiLayer
from mpfit.mpfit3 import mpfit


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
    psfdict = {"mGL":3, "3DG":1, "GL":0}
    psf_nparasdict = {"mGL":3, "3DG":2, "GL":2}


    def __init__(self, psfmodel="mGL", zoffset=None, channels=[1]):
        super(zscanFitter, self).__init__()
        if psfmodel in self.psfdict.keys():
            self._psfmodel = psfmodel     # psfmodel
        else:
            raise ValueError("{} is not a correct PSF model".format(psfmodel))

        self._zoffset = zoffset

        if not isinstance(channels, list):
            raise TypeError("The type of channels is a list.")
        self._channels = channels
        # psf can differ from channel to channel due to color abberation.
        self._psf, self._psffixed = {}, {}
        for cha in channels:
            self._psf.setdefault(cha, []) # for psf paras on each channel
            self._psffixed.setdefault(cha, []) # fixed condition for psf
        # geoinfo is given as a list of layers
        self._geoinfo = [] # for sample geometric models
        self._geoinfofixed = [] # fixed for sample geometric models
        # spillover paras for multi-channel zscan FFS
        self._spillover = {}

    def setkz(self, z, kz=[]):
        # TODO
        # self.plot()
        pass

    def _checkfitinputs(self, z, kz, errkz):
        """Check if the fit inputs are correct.


        """
        if len(kz) != len(self.channels):
            raise ValueError("The number of kzs does not match to the number \
                            of channels")

        if (errkz is not None):
            if not isinstance(errkz, list):
                raise TypeError("The type of errkz is not list")
            elif len(errkz) != len(self.channels):
                raise ValueError("The number of kzs does not match to \
                                the number of channels")

        for cha in self.channels:
            if self._isPSFEmpty(cha):
                raise ValueError("The psf paras for the channel {} is not \
                                given.".format(cha))
        # check the zoffset
        if self._isZoffsetEmpty():
            tzoffset = input("Set up zoffset in the range between \
                    {%.2f} and {%.2f} : ".format(z.min(), z.max()))
            self._zoffset = float(tzoffset)
        elif not (z.min() <= self._zoffset <= z.max()):
            print("{} is out of the range between \
                    {%.2f} and {%.2f} : ".format(z.min(), z.max()))
            tzoffset = input("Set up zoffset in the range between \
                    {%.2f} and {%.2f} : ".format(z.min(), z.max()))
            self._zoffset = float(tzoffset)

        # check the geo info
        if self._isGeoinfoEmpty():
            raise ValueError("No geoinfo is available.")

        # check the spillover
        if len(self.channels) > 1 and self._isSpillOverEmpty():
            raise ValueError("No spillover paras is available.")


    @staticmethod
    def getErrkz(kz):
        """Calculate the error for kz.
        sqrt(phton count) is equal to the error of photon count due to the shot
        noise.
        """
        if not isinstance(kz, list):
            raise TypeError("The input kz is not list")
        res = [np.sqrt(x) for x in kz]
        for x in res:
            index = (x == 0)
            x[index] += 1.
        return res

    def fit(self, z, kz=[], errkz=None):
        """fit zscan intensity profiles of multi-channel with a single geometric
        model

        Questions: How to apply the spillover parameters on each each geometric
        layer. especially for the background counts when the fluorescent
        intensity is very low.
        """
        self._checkfitinputs(z, kz, errkz)

        if errkz is None:
            errkz = self.getErrkz(kz)

        # TODO
        x = z
        y = np.array(kz).flatten()
        yerr = np.array(errkz).flatten()
        paras, fixed, fitinfo = self._generateparas()
        # see docstring in _generateparas to know what paras, fitinfo, and fixed
        # are.

        parinfo = [{'value':v, 'fixed':f, 'limited':[1,0], 'limits':[0.,0.]}
         								for v, f in zip(paras, fixed)]

        # TODO consider extra features in myfunct
        def myfunct(p, fjac=None, x=None, y=None, err=None, info=None):
            model =  self.kzMultiLayerFCT(x, p, info=info)
            status = 0
            return [status, (y-model)/err]

        fa = {"x":x, "y":y, "err":yerr, "info":fitinfo}
        res = mpfit(myfunct, paras, functkw=fa, parinfo=parinfo, maxiter=300, quiet=1)
        yfit = self.kzMultiLayerFCT(x, res.params, info=fitinfo)
        return res, yfit

    def _checkpsfmodel(self, psfparas, fixed):
        if self._psfmodel == "mGL":
            if psfparas != []:
                assert (len(psfparas) ==  3), "The number of paras should be 3."
        elif self._psfmodel == "GL":
            if psfparas != []:
                assert (len(psfparas) ==  2), "The number of paras should be 2."
        elif self._psfmodel == "3DG":
            if psfparas != []:
                assert (len(psfparas) ==  2), "The number of paras should be 2."

        if fixed != []:
            assert (len(psfparas) ==  len(fixed)), \
                "The number of paras does not match to the number of fixed."

    def setPSF(self, channel=1, psfparas=[], fixed=[]):
        """set PSF model's paras and fixed conditions

        psfparas = [zR, y, w0] for mGL
        psfparas = [zR, w0] for GL
        psfparas = [zR, w0] for 3DG

        """
        if channel not in self._channels:
            raise ValueError("The channel {} is not in channels".format(channel))

        self._checkpsfmodel(psfparas, fixed)
        if psfparas != [] and fixed != []:
            self._psf[channel] = psfparas
            self._psffixed[channel] = fixed
        elif psfparas != [] and fixed == []:
            self._psf[channel] = psfparas
            self._psffixed[channel] = [0]*len(psfparas)
        elif psfparas == [] and fixed == []:
            self._setPSFparas(channel)
        else:
            raise ValueError("psfpara is not available.")

    def _setPSFparas(self, channel):
        """set PSF paras by hands"""
        print("set PSF paras for {} psf model".format(self._psfmodel))
        self._psf.setdefault(channel, []) # for psf paras on each channel
        self._psffixed.setdefault(channel, []) # fixed condition for psf

        if self._psfmodel == "mGL":
            psfpara_names = ['zR', 'y', 'w0']
        elif self._psfmodel in ['GL', '3DG']:
            psfpara_names = ['zR', 'w0']

        for i in psfpara_names:
            h1 = input("PSF para for {} ? ".format(i))
            h2 = input("fix {} (1 or 0) ? ".format(i))
            self._psf[channel].append(float(h1))
            self._psffixed[channel].append(float(h1))

    def getPSF(self, channel=1):
        """get the whole information about the PSF"""
        return self._psfmodel, self._psf[channel], self._psffixed[channel]

    def addLayer(self, geomodel):
        """geomodel is either integer or string in geodict.
        paras and fixed are manually set up.

        geomodel >> a key in geodict or a value in geodict
        """
        self.setLayer(geomodel)
        self._setGeoParas()
        return

    def _setGeoParas(self):
        """set geo-paras and fixed by hands"""
        print("set paras for {} model".format(self._geoinfo[-1]['geo']))
        para_names, fpara_names = self._paranames()
        for i, j in zip(para_names, fpara_names):
            h1 = input("para for {} ? ".format(i))
            h2 = input("fix {} (1 or 0) ? ".format(i))
            self._geoinfo[-1][i] = float(h1)
            self._geoinfofixed[-1][j] = int(h2)
        return

    def _checkgeomodel(self, geomodel):
        if isinstance(geomodel, int):
            for key, value in self.geodict.items():
                if value == geomodel:
                    return key
            else:
                raise ValueError("geomodel is not available.")
        elif isinstance(geomodel, str):
            if geomodel in self.geodict:
                return geomodel
            else:
                raise ValueError("geomodel is not available.")
        else:
            raise TypeError("The type of geomodel is incorrect.")

    def _paranames(self):
        if len(self.channels) == 1:
            return ['k1', 'LR'], ['fk1', 'fLR']
        elif len(self.channels) == 2:
            return ['k1', 'k2', 'LR'], ['fk1', 'fk2', 'fLR']
        elif len(self.channels) == 3:
            return ['k1', 'k2', 'k3', 'LR'], ['fk1', 'fk2', 'fk3', 'fLR']
        else:
            raise ValueError("The analysis is not available.")

    def setLayer(self, geomodel, paras=[], fixed=[], layer_index=None):
        """set geometric model on each layer with geo-paras and fixed conditions

        geomodel : a key in geodict or a value in geodict

        paras: list
        for single channel [1]        : [k1, LR]
        for dual channels [1, 2]      : [k1, k2, LR]
        for triple channels [1, 2, 3] : [k1, k2, k3, LR]
        where k1, k2, k3 >> counts per bin, LR >> length or radius

        fixed: list of 0 or 1 (0: free, 1: fixed in fitting)
        For a given fixed, len(fixed) == len(paras)

        layer_index : None or a non-negative integer
        None : a single geometric layer is added at the end of current
        geometric models
        a non-negative integer : a geometric layer of the layer_index is reset
        by given paras and fixed

        """
        geo = self._checkgeomodel(geomodel)

        para_names, fpara_names = self._paranames()
        if paras == [] and fixed == []:
            paras = [0.]*len(para_names)
            fixed = [0]*len(fpara_names)
        elif paras != [] and fixed == []:
            fixed = [0]*len(fpara_names)
        elif paras == [] and fixed != []:
            paras = [0.]*len(para_names)
        else:
            assert len(paras) == len(para_names), "The number elements in paras\
                    should be equal to {}".format(len(self.channels) + 1)
            assert len(fixed) == len(fpara_names), "The number elements in fixed\
                    should be equal to {}".format(len(self.channels) + 1)

        if len(paras) != len(fixed):
            raise ValueError("The number of elements in paras does not \
                        match to the number of elements in fixed.")
        # assign paras and fixed to _geoinfo and _geoinfofixed
        if layer_index == None:
            self._geoinfo.append(dict([(k, v) for k, v in zip(para_names, paras)]))
            self._geoinfofixed.append(dict([(k, v) for k, v in zip(fpara_names, fixed)]))
            self._geoinfo[-1]['geo'] = geo
            self._geoinfofixed[-1]['geo'] = geo
        elif (isinstance(layer_index, int) and
                layer_index < len(self._geoinfo)):
            if self._geoinfo[layer_index]['geo'] == geo:
                for k, v in zip(para_names, paras):
                    self._geoinfo[layer_index][k] = v
                for k, v in zip(fpara_names, fixed):
                    self._geoinfofixed[layer_index][k] = v
            else:
                raise ValueError("model does not match.")
        else:
            raise ValueError("layer_index is out of the range")
        return

    def removeLayer(self, layer_index=None):
        # remove a geometric layer
        if self._geoinfo != []:
            try:
                if isinstance(layer_index, int):
                    self._geoinfo.pop(layer_index)
                    self._geoinfofixed.pop(layer_index)
                elif layer_index == None:
                    self._geoinfo.pop()
                    self._geoinfofixed.pop()
            except:
                raise ValueError("layer_index is out of the allowed range")
        else:
            print("geoinfo is empty.")
        return

    def setSpillover(self, paras):
        """a method for setting up spillover parameters
        """
        if len(self.channels) == 1:
            return
        elif len(self.channels) == 2:
            spillover = ['f12']
        elif len(self.channels) == 3:
            spillover = ['f12', 'f13', 'f23']
        else:
            raise ValueError("The analysis is not available")
        assert len(paras) == len(spillover), "The number of elments in paras\
            should be equal to {}".format(len(spillover))
        for k, v in zip(spillover, paras):
            self._spillover[k] = v
        return



    def _generateparas(self):
        """
        paras = [psfparas for channel 1, (psfparas for channel 2, psfparas for channel 3),
                offset, layer[0]_paras, layer[1]_paras, ......]
        fixed =
        fitinfo = {"nch":#, "psfmodel":#, "n_psfparas":#, "geo":[]}
        """
        for cha in self._channels:
            if self._isPSFEmpty(cha):
                raise ValueError("psf paras are not available.")
        if self._isGeoinfoEmpty():
            raise ValueError("geoinfo paras are not available.")
        if self._isZoffsetEmpty():
            raise ValueError("zoffset para is not available.")

        paras, fixed, fitinfo = [], [], {}
        nch = 0
        for channel in sorted(self._psf.keys()):
            paras += self._psf[channel]
            fixed += self._psffixed[channel]
            nch += 1
        fitinfo["nch"] = nch
        fitinfo["psfmodel"] = self.psfdict[self._psfmodel]
        fitinfo["n_psfparas"] = self.psf_nparasdict[self._psfmodel]

        fitinfo["geo"] = []
        paras += [self._zoffset]
        fixed += [0]
        para_names, fpara_names = self._paranames()
        for i, j in zip(self._geoinfo, self._geoinfofixed):
            temp = [i[k] for k in para_names]
            ftemp = [j[k] for k in fpara_names]
            paras.extend(temp)
            fixed.extend(ftemp)
            fitinfo["geo"] += [self.geodict[i["geo"]]]

        fitinfo["spillover"] = self._spillover

        return np.array(paras).flatten(), fixed, fitinfo

    # def kzfct(self, channel=1):
    #     self._geoinfo
    #     zscanMultiLayer
    #     a, b = self._generateparas()
    #     kz_fct = self.kzMultiLayerFCT(self.z, a, info=b)

    @staticmethod
    def kzMultiLayerFCT(z, paras, info=None):
        """zscan multilayer function for the fit.

        For given paras and info,
            parasPSf, zoffset, model are reconstituted for zscanMultiLayer.

        """
        result = np.zeros(z.size)  # z.size / nch

        nch = info["nch"]
        n_psfp0 = info["n_psfparas"]

        psfparas = []
        for i in range(nch):
            psfparas.append(paras[0 + i*n_psfp0: n_psfp0 + i*n_psfp0])
        zoff = paras[n_psfp0*nch]

        nparas = {1:2, 2:3, 3:4}   # nparas[nch] == len(self._paranames()[0])
        if nch == 2:
            spo = ['f12']
        elif nch == 3:
            spo = ['f12', 'f13', 'f23']

        geomodels = [[] for x in range(nch)]
        for x in range(len(info["geo"])):
            temp = paras[n_psfp0*nch + 1 + nparas[nch]*x
                            :n_psfp0*nch + nparas[nch] + 1 + nparas[nch]*x]
            geomodels[0].append({"geo":info["geo"][x], "k":temp[0], "LR":temp[nparas[nch]-1]})
            if nch == 2:
                geomodels[1].append({"geo":info["geo"][x], "k":temp[0], "LR":temp[nparas[nch]-1]})
            elif nch == 3:
                geomodels[1].append({"geo":info["geo"][x], "k":temp[1], "LR":temp[nparas[nch]-1]})
                geomodels[2].append({"geo":info["geo"][x], "k":temp[2], "LR":temp[nparas[nch]-1]})

        if nch == 1:
            return zscanMultiLayer(z, zoff, psfparas[0], model=geomodels[0], psfmodel=info["psfmodel"])
        else:
            # TODO zscan profiles for multiple channel
            # TODO take into account the spillover.
            res = []
            for i in range(nch):
                t = zscanMultiLayer(z, zoff, psfparas[i], model=geomodels[i], psfmodel=info["psfmodel"])
                if i == 0:
                    res.append(t)
                elif i == 1:
                    a = info["spillover"][spo[0]]
                    t += a*res[0]
                    res.append(t)
                elif i == 2:
                    a = info["spillover"][spo[1]]
                    b = info["spillover"][spo[2]]
                    res.append()
            temp = np.hstack(res)
            # print(temp)
            return temp
            # return np.concatenate(res, axis=1)


    def _isPSFEmpty(self, channel):
        return self._psf[channel] == []

    def _isGeoinfoEmpty(self):
        return self._geoinfo == []

    def _isZoffsetEmpty(self):
        return (self._zoffset is None)

    def _isSpillOverEmpty(self):
        return self._spillover == {}


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

    @property
    def zoffset(self):
        return self._zoffset
    @zoffset.setter
    def zoffset(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("The type of zoffset value is either int or float.")
        else:
            self._zoffset = value

    @property
    def channels(self):
        return self._channels
    @channels.setter
    def channels(self, values):
        if not isinstance(values, list):
            raise TypeError("The type of channels is list.")
        else:
            self._channels = values

    @property
    def spillover(self):
        return self._spillover

    @classmethod
    def printgeodict(cls):
        temp = sorted(cls.geodict.items(), key= (lambda x:x[1]))
        for key, value in temp:
            print(key, value)

def main():
    from zscanTransformer import zscanTransformer as zscan
    from readFFSfromFLEX import readFFSfromFLEX as ffs
    from matplotlib import pyplot as plt

    data = ffs(["zscan_slab_egfp.dat"], [1, 2], 20000)
    temp_zscan = zscan(channels=[2], slice_zscans = True)
    res = temp_zscan.transform(data)

    zscanfit = zscanFitter(psfmodel="mGL", zoffset=13., channels=[1])
    zscanfit.setPSF(channel=1, psfparas=[1., 2., 0.45], fixed=[0, 0, 0])
    # print("psf :", zscanfit.getPSF())
    #zscanfit.addLayer("DOWN")
    zscanfit.setLayer("DOWN", [1., 0.], [0, 1])
    zscanfit.setLayer(0, [1300., 1.], [0, 0])
    zscanfit.setLayer("UP", [1., 0.], [0, 1])
    #zscanfit.setLayer("UP", [5., 5.], [0, 0], layer_index=1)

    zz, yfit = zscanfit.fit(res[0], [res[2][0]])
    print(zz.params)
    # for x in res[2]:
    #     plt.plot(res[0], x)
    plt.plot(res[0], res[2][0])
    plt.plot(res[0], yfit, 'r')
    plt.xlabel('z (um)')
    plt.ylabel('counts per {} bins'.format(temp_zscan.nbins))
    plt.show()


    temp_zscan = zscan(channels=[1, 2], slice_zscans = True)
    res = temp_zscan.transform(data)
    # for x, y in zip(res[2], res[1]):
    #     plt.plot(res[0], x)
    #     plt.plot(res[0], y)
    plt.plot(res[0], res[2][0])
    plt.plot(res[0], res[1][0])
    plt.show()

    zscanfit_dual = zscanFitter(psfmodel="mGL", zoffset=13., channels=[1,2])
    zscanfit_dual.setPSF(channel=1, psfparas=[1., 2., 0.45], fixed=[0, 0, 0])
    zscanfit_dual.setPSF(channel=2, psfparas=[1., 2., 0.45], fixed=[0, 0, 0])
    zscanfit_dual.setSpillover([1./8.])
    zscanfit_dual.setLayer("DOWN", [1., 1., 0.], [0, 0, 1])
    zscanfit_dual.setLayer(0, [1300., 400., 1.], [0, 0, 0])
    zscanfit_dual.setLayer("UP", [1., 1., 0.], [0, 0, 1])

    zz, yfit = zscanfit_dual.fit(res[0], [res[2][0], res[1][0]])
    print(zz.params)
    v = yfit.reshape((2, res[0].size))
    print(v)
    print(v.shape)
    plt.plot(res[0], res[2][0])
    plt.plot(res[0], res[1][0])
    plt.plot(res[0], v[0], 'r')
    plt.plot(res[0], v[1], 'r')
    plt.show()


if __name__=="__main__":
    main()
