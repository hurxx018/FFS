from scipy.special import gamma, hyp2f1, erf
from scipy.integrate import quad
import numpy as np

def erfapprox(x):
    "approximation of erf function"
    a = 0.14001229
    x2 = np.power(x, 2)
    q = (4./np.pi + a*x2) / (1. + a*x2)
    sgnx = (x > 0)*2 - 1
    return sgnx * np.sqrt(1. - np.exp(-x2 * q))

def volPSFGLspCylinderGeoIntegrand(zeta, dz, r2, b, s):
    #; p = {dz:dz, R2:R2, b:b, s:s}
    zeta2p1 = 1. + np.power(zeta, 2)
    return erf(2.*b*np.sqrt(r2-np.power(dz + zeta, 2))/np.sqrt(zeta2p1)) \
               / np.power(zeta2p1, s)
    # return erfapprox(2.*b*np.sqrt(r2- np.power(dz + zeta, 2))/np.sqrt(zeta2p1)) \
    #         / np.power(zeta2p1, s)

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
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2]) #beam waist ratio b=zR/w0
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
    b   = np.float64(parasPSF[0])/np.float64(parasPSF[2]) #beam waist ratio b=zR/w0
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


def integralPSFglsp(x, y):
    # calculate the integration of modified GL psf
    tx = np.float64(x)
    return tx*hyp2f1(0.5, y, 1.5, -np.power(tx, 2))

def volPSFGLsp(paras, z1, z2, geo=None, absolute=False):
    """
    calculate normalized volume of GL PSF restricted to the z-interval [z1,z2]
    with PSF centered at (z2 + z1)/2

    if KEYWORD absolute is set, the pgm calculates the regular volume
    KEYWORD GEO defines the gemoetry of the sample

    INPUT
    paras  [zR, y, w0] = (axial waist, index, radial waist)
    #paras  [zR, y]     = (axial waist, index)                 if absolute=True
    z1     lower z-limit
    z2     upper z-limit

    KEYWORDS
    geo =   0 = Slab          (Default)
    geo =   1 = SemiInfiniteUp
    geo =   2 = SemiInfiniteDown
    geo =   3 = Delta function
    geo =  -1 = Infinite

    absolute = 0  (Default)
    absolute = 1

    RETURN
    normalized volume of PSF

    Call sequence

    VolPSFGLsp, [zR,y    ], z1, z2 {,Geo = 0}  --> Vol[  z1,   z2] / VolInfinity
    VolPSFGLsp, [zR,y    ], z1, z2  ,Geo = 1   --> Vol[  z1, +Inf] / VolInfinity
    VolPSFGLsp, [zR,y    ], z1, z2  ,Geo = 2   --> Vol[-Inf,   z2] / VolInfinity
    VolPSFGLsp, [zR,y    ], z1      ,Geo = 3   --> Area[z1] / AreaMax
    VolPSFGLsp, [zR,y    ]          ,Geo =-1   --> 1

    VolPSFGLsp, [zR,y, w0], z1, z2 {,Geo= 0} ,absolute=1 --> Vol[  z1,   z2]
    VolPSFGLsp, [zR,y, w0], z1, z2  ,Geo= 1  ,absolute=1 --> Vol[  z1, +Inf]
    VolPSFGLsp, [zR,y, w0], z1, z2  ,Geo= 2  ,absolute=1 --> Vol[-Inf,   z2]
    VolPSFGLsp, [zR,y, w0], z1      ,Geo= 3  ,absolute=1 --> Area[z1]
    VolPSFGLsp, [zR,y, w0]          ,Geo=-1  ,absolute=1 --> Vol[-Inf, +Inf] = VolInfinity
    """

    if geo == None:
        geometry = 0
    else:
        geometry = geo

    zR  = np.float64(paras[0])
    y   = np.float64(paras[1])
    VolInf = np.float64(1.)
    fac = gamma(y)/gamma(y - 0.5)/np.sqrt(np.pi)

    if absolute:
        w0  = np.float64(paras[2])
        w02 = np.power(w0, 2)
        VolInf = np.pi/4. *w02*zR / fac
        if geometry == 3:
            VolInf = np.pi/4. * w02   #infinite area

    if geometry == 0:
        vol = (integralPSFglsp(z2/zR, y) - integralPSFglsp(z1/zR, y)) * fac
    elif geometry == 1:
        vol = 0.5 - integralPSFglsp(z1/zR, y) * fac
    elif geometry == 2:
        vol = 0.5 + integralPSFglsp(z2/zR, y) * fac
    elif geometry == 3:
        vol = 1./np.power(1. + np.power(z1/zR, 2.), y)              # area
    elif geometry == -1:
        vol = 1.
    else:
        raise ValueError("The geometry is not available.")

    return vol * VolInf

def zscanPSFGLsp(z, parasPSF, parasGEO, geo=None, absolute=False):
    """
    calculate the normalized PSF-volume as a function of z
            for a given sample geometry

    KEYWORDS:
    geo   defines geometry of sample
    absolute  if set, the regular volume is calculated

    INPUT:
    z           z values (vector)

    parasPSF   [zR, y]     = (axial waist, index)           if Geo = {0,1,2,3} AND absolute=0
    parasPSF   [zR, y, w0] = (axial waist, index, radial waist) if Geo = {0,1,2,3} AND absolute=1
    parasPSF   [zR, y, b ] = (axial waist, index, b = zR/w0)    if Geo = {4,5}

    parasGEO   [zoff, L] = (z-offset to center of slab, Length of slab) if GEO =   0 = Slab
    parasGEO   [zoff]    = (z-offset to beginning of SemiInfinite Slab) if Geo =   1 = SemiInfiniteUp
    parasGEO   [zoff]    = (z-offset to end of SemiInfinite Slab)       if Geo =   2 = SemiInfiniteDown
    parasGEO   [zoff]    = (z-offset to delta layer)                    if Geo =   3 = Delta
    parasGEO   [zoff, R] = (z-offset to center, radius)                 if Geo =   4 = Sphere
    parasGEO   [zoff, R] = (z-offset to center, radius)                 if Geo =   5 = Cylinder


    geo:
       0 = slab          (default)
       1 = semiinfinite up
       2 = semiinfinite down
       3 = delta
       4 = sphere
       5 = cylinder

    absolute:
       False = normalized volume (default)
       True = absolute volume

    RETURN
    volume of PSF as a function of z

    Call sequence
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff,L]=parasGEO, {,Geo = 0}  --> normalized volume @ zexp for Slab
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff, ]=parasGEO,  ,Geo = 1   --> normalized volume @ zexp for SemiInfiniteUp
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff, ]=parasGEO,  ,Geo = 2   --> normalized volume @ zexp for SemiInfiniteDown
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff, ]=parasGEO,  ,Geo = 3   --> normalized volume @ zexp for Delta
    zscanPSFGLsp(z, [zR,y,b]=parasPSF, [zoff,R]=parasGEO,  ,Geo = 4   --> normalized volume @ zexp for Spherical Geometry
    zscanPSFGLsp(z, [zR,y,b]=parasPSF, [zoff,R]=parasGEO,  ,Geo = 5   --> normalized volume @ zexp for Cylindrical Geometry

    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff,L]=parasGEO, {,geo = 0}, absolute = 1  --> absolute volume @ zexp for Slab
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff  ]=parasGEO,  ,geo = 1 , absolute = 1  --> absolute volume @ zexp for SemiInfiniteUp
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff  ]=parasGEO,  ,geo = 2 , absolute = 1  --> absolute volume @ zexp for SemiInfiniteDown
    zscanPSFGLsp(z, [zR,y  ]=parasPSF, [zoff  ]=parasGEO,  ,geo = 3 , absolute = 1  --> absolute volume @ zexp for Delta
    zscanPSFGLsp(z, [zR,y,b]=parasPSF, [zoff,R]=parasGEO,  ,geo = 4 , absolute = 1  --> absolute volume @ zexp for Spherical Geometry
    zscanPSFGLspv(z, [zR,y,b]=parasPSF, [zoff,R]=parasGEO,  ,geo = 5 , absolute = 1  --> absolute volume @ zexp for Cylindrical Geometry

    """

    if geo == None:
        geometry = 0
    else:
        geometry = geo

    zR   = np.float64(parasPSF[0])
    y    = np.float64(parasPSF[1])
    w0   = np.float64(parasPSF[2])
    zoff = np.float64(parasGEO[0])
    zPSF = (z - zoff).astype(np.float64)

    if geometry == 0:
        L   = parasGEO[1]
        vol = volPSFGLsp(parasPSF, -L/2. - zPSF, L/2. - zPSF, absolute = absolute)
    elif geometry == 1:
        vol = volPSFGLsp(parasPSF, - zPSF, 0, geo = geometry, absolute = absolute)
    elif geometry == 2:
        vol = volPSFGLsp(parasPSF, 0, -zPSF, geo = geometry, absolute = absolute)
    elif geometry == 3:
        vol = volPSFGLsp(parasPSF, zPSF, None, geo = geometry, absolute = absolute)
    elif geometry == 4:
        R = parasGEO[1]
        vol = volPSFGLspSphereGeo(parasPSF, zPSF, R, absolute = absolute)
    elif geometry == 5:
        R = parasGEO[1]
        vol = volPSFGLspCylinderGeo(parasPSF, zPSF, R, absolute = absolute)
    else:
        raise ValueError("geometry is not available.")

    return vol
