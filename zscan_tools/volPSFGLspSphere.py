from scipy.special import gamma
from scipy.integrate import quad
import numpy as np

def volPSFGLspSphereGeoIntegrand(zeta, dz, r2, b2, s):
    zeta2p1 = 1. + np.power(zeta, 2)
    return (1. - np.exp(-4.*b2*(r2 - np.power(zeta + dz, 2))/zeta2p1)) \
                    / np.power(zeta2p1, s)


def volPSFGLspSphereGeo(parasPSF, dz, R, absolute=False):
    """
    calculate normalized volume of modified GL PSF for a spherical sample geometry with beam along z-axis
    (with z-axis passing through center of sphere)
    Lab coordinate system origin at center of sphere (dz=0)
    with PSF centered at 0

    if KEYWORD absolute is set, the pgm calculates the regular volume

    INPUT
    parasPSF  [zR, y, w0] = (axial waist, index, radial waist w0)
    dz     center of PSF along z-axis
    R      Radius of sphere

    absolute = False  (Default)
    Absolute = True

    RETURN
    normalized volume of PSF

    Call sequences
    volPSFGLspSphereGeo, [zR,y,w0], dz, R            --> Vol[  dz,   R] / VolInfinity
    volPSFGLspSphereGeo, [zR,y,w0], dz, R, absolute=True    --> Vol[  dz,   R]
    """
    zR  = np.float64(parasPSF[0])
    y   = np.float64(parasPSF[1])
    #beam waist ratio b=zR/w0
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2])
    fac = gamma(y)/gamma(y-0.5)/np.sqrt(np.pi)

    if absolute:
        w0  = np.float64(parasPSF[2])
        w02 = np.power(w0, 2.)
        fac = np.pi/4. *w02*zR

    vol = np.zeros_like(dz, dtype=np.float64)

    Rt  = R/zR
    dzt = dz/zR
    z1 = -dzt - Rt
    z2 = -dzt + Rt
    Rt2 = np.power(Rt, 2.)
    b2  = np.power(b, 2.)

    if dzt.size == 1:
        vol = quad(volPSFGLspSphereGeoIntegrand, z1, z2,
                        args=(dzt, Rt2, b2, y))[0]
    else:
        for i, val in enumerate(dzt):
            vol[i] = quad(volPSFGLspSphereGeoIntegrand, z1[i], z2[i],
                        args=(val, Rt2, b2, y))[0]
    return vol * fac

def volPSFsqrGLspSphereGeo(parasPSF, dz, R, absolute=False):
    """
    calculate normalized squared volume of modified GL PSF for a spherical sample geometry with beam along z-axis
    (with z-axis passing through center of sphere)
    Lab coordinate system origin at center of sphere (dz=0)
    with PSF centered at 0

    if KEYWORD absolute is set, the pgm calculates the regular volume

    INPUT
    parasPSF  [zR, y, w0] = (axial waist, index, radial waist w0)
    dz     center of PSF along z-axis
    R      Radius of sphere

    absolute = False  (Default)
    Absolute = True

    RETURN
    normalized volume of PSF

    Call sequences
    volPSFsqrGLspSphereGeo, [zR,y,w0], dz, R            --> Vol[  dz,   R] / VolInfinity
    volPSFsqrGLspSphereGeo, [zR,y,w0], dz, R, absolute=True    --> Vol[  dz,   R]
    """
    zR  = np.float64(parasPSF[0])
    y   = np.float64(parasPSF[1])
    #beam waist ratio b=zR/w0
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2])
    fac = 0.5*gamma(y)/gamma(y-0.5)/np.sqrt(np.pi)

    if absolute:
        w0  = np.float64(parasPSF[2])
        w02 = np.power(w0, 2.)
        fac = 0.5*np.pi/4. *w02*zR

    vol = np.zeros_like(dz, dtype=np.float64)

    Rt  = R/zR
    dzt = dz/zR
    z1 = -dzt - Rt
    z2 = -dzt + Rt
    Rt2 = np.power(Rt, 2.)
    b2  = np.power(b, 2.)

    if dzt.size == 1:
        vol = quad(volPSFGLspSphereGeoIntegrand, z1, z2,
                        args=(dzt, Rt2, 2*b2, 2*y+1))[0]
    else:
        for i, val in enumerate(dzt):
            vol[i] = quad(volPSFGLspSphereGeoIntegrand, z1[i], z2[i],
                        args=(val, Rt2, 2*b2, 2*y+1))[0]
    return vol * fac
