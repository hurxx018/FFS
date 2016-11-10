from scipy.special import gamma, erf
from scipy.integrate import quad
import numpy as np

def erfapprox(x):
    "approximation of erf function"
    a = 0.14001229
    x2 = np.power(x, 2)
    q = (4./np.pi + a*x2) / (1. + a*x2)
    sgnx = (x > 0)* - 1
    return sgnx * np.sqrt(1. - np.exp(-x2 * q))

def volPSFGLspCylinderGeoIntegrand(zeta, dz, r2, b, s):
    #; p = {dz:dz, R2:R2, b:b, s:s}
    zeta2p1 = 1. + np.power(zeta, 2)
    #return, erf(2.*b*np.sqrt(r2-np.power(dz + zeta, 2))/np.sqrt(zeta2p1)) \
    #            / np.power(zeta2p1, s)
    return erfapprox(2.*b*np.sqrt(r2- np.power(dz + zeta, 2))/np.sqrt(zeta2p1)) \
            / np.power(zeta2p1, s)

def volPSFGLspCylinderGeo(parasPSF, dz, R, absolute=False):
    """
    calculate normalized volume of modified GL PSF for a cylindrical sample
    geometry with beam along z-axis
    assumes very long cylinder with long axis lying in the xy-plane.
    Z-axis passes through the center of the cylinder's crosssection (disk)
    with Lab coordinate origin at center of disk (dz=0) with PSF centered at dz

    if KEYWORD absolute is set to be True, the pgm calculates the regular volume

    INPUT
    parasPSF  [zR, y, w0] = (axial waist, index, radial beam waist)
    dz     center of PSF along z-axis
    R      Radius of cylinder

    absolute = False  (Default)
    absolute = True

    RETURN
    normalized volume of PSF

    Call sequence
    volPSFGLspCylinderGeo, [zR,y,w0], dz, R               --> Vol[  dz,   R] / VolInfinity
    volPSFGLspCylinderGeo, [zR,y,w0], dz, R,Absolute=1    --> Vol[  dz,   R]
    """
    zR  = np.float64(parasPSF[0])
    y   = np.float64(parasPSF[1])
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2]) beam waist ratio b=zR/w0
    fac = gamma(y)/gamma(y - 0.5)/np.sqrt(np.pi)

    if absolute:
        w0  = np.float64(parasPSF[2])
        w02 = np.power(w0, 2)
        fac = np.pi/4. *w02*zR

    vol = np.zeros_like(dz, dtype=np.float64)

    Rt  = R/zR
    dzt = dz/zR
    z1 = -dzt - Rt
    z2 = -dzt + Rt
    Rt2 = np.power(Rt, 2)

    if dzt.size == 1:
        vol = quad(volPSFGLspCylinderGeoIntegrand, z1, z2,
                        args=(dzt, Rt2, b, y))[0]
    else:
        for i, val in enumerate(dzt):
            vol[i] = quad(volPSFGLspCylinderGeoIntegrand, z1[i], z2[i],
                        args=(val, Rt2, b, y))[0]
    return vol * fac

def volPSFsqrGLspCylinderGeo(parasPSF, dz, R, absolute=False):
    """
    calculate normalized volume of squared modified GL PSF
    for a cylindrical sample geometry with beam along z-axis
    assumes very long cylinder with long axis lying in the xy-plane.
    Z-axis passes through the center of the cylinder's crosssection (disk)
    with Lab coordinate origin at center of disk (dz=0) with PSF centered at dz

    if KEYWORD absolute is set to be True, the pgm calculates the regular volume

    INPUT
    parasPSF  [zR, y, w0] = (axial waist, index, radial beam waist)
    dz     center of PSF along z-axis
    R      Radius of cylinder

    absolute = False  (Default)
    absolute = True

    RETURN
    normalized volume of PSF

    Call sequence
    volPSFsqrGLspCylinderGeo, [zR, y, w0], dz, R --> Vol[  dz,   R] / VolInfinity
    volPSFsqrGLspCylinderGeo, [zR, y, w0], dz, R,Absolute=1   --> Vol[  dz,   R]
    """
    zR  = np.float64(parasPSF[0])
    y   = np.float64(parasPSF[1])
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2]) beam waist ratio b=zR/w0
    fac = 0.5*gamma(y)/gamma(y - 0.5)/np.sqrt(np.pi)

    if absolute:
        w0  = np.float64(parasPSF[2])
        w02 = np.power(w0, 2)
        fac = 0.5*np.pi/4. *w02*zR

    vol = np.zeros_like(dz, dtype=np.float64)

    Rt  = R/zR
    dzt = dz/zR
    z1 = -dzt - Rt
    z2 = -dzt + Rt
    Rt2 = np.power(Rt, 2)

    if dzt.size == 1:
        vol = quad(volPSFGLspCylinderGeoIntegrand, z1, z2,
                        args=(dzt, Rt2, np.sqrt(2.)*b, 2.*y+1.))[0]
    else:
        for i, val in enumerate(dzt):
            vol[i] = quad(volPSFGLspCylinderGeoIntegrand, z1[i], z2[i],
                        args=(val, Rt2, np.sqrt(2.)*b, 2.*y+1.))[0]
    return vol * fac
