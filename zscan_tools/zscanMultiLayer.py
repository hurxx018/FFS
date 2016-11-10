import numpy as np


def zscanPSF(z, parasPSF, parasGEO, geo=None, psfmodel=None, order=None):
    """
    calculate the normalized PSF volume of a PSF model (mGL, 3DG, GL)
        scanned through a sample

    z          z values (vector)
    parasPSF   [zR,y]    = (axial waist, y           for GEO = {0,1,2,3})
    parasPSF   [zR,y,b]  = (axial waist, y, b=zR/w0  for GEO = {4, 5})

    parasGEO   [zoff, L] = (z-offset to center of slab, Length of slab) if GEO =   0 = Slab
    parasGEO   [zoff   ] = (z-offset to beginning of SemiInfinite Slab) if Geo =   1 = SemiInfiniteUp
    parasGEO   [zoff   ] = (z-offset to end of SemiInfinite Slab)       if Geo =   2 = SemiInfiniteDown
    parasGEO   [zoff   ] = (z-offset to delta layer)                    if Geo =   3 = Delta
    parasGEO   [zoff, R] = (z-offset from center of sphere,   radius)   if Geo =   4 = Sphere
    parasGEO   [zoff, R] = (z-offset from center of cylinder, radius)   if Geo =   5 = Cylinder
    parasGEO   [zoff, L] = (z-offset to center of slab, Length of slab) if Geo =  10 = NotSlab
    parasGEO   [zoff   ] = (z-offset to delta layer)                    if Geo =  13 = NotDelta
    parasGEO   [zoff, R] = (z-offset from center of sphere,   radius)   if Geo =  14 = NotSphere
    parasGEO   [zoff, R] = (z-offset from center of cylinder, radius)   if Geo =  15 = NotCylinder
    parasGEO   [zoff, R] = (z-offset from center of sphere,   radius)   if Geo =  24 = NotSphereUp
    parasGEO   [zoff, R] = (z-offset from center of cylinder, radius)   if Geo =  25 = NotCylinderUp

    KEYWORDS INPUT
    GEO:
    0 = slab          (default)
    1 = semiinfinite up
    2 = semiinfinite down
    3 = delta
    4 = sphere
    5 = cylinder
    10 = NOTslab
    13 = NOTdelta
    14 = NOTsphere
    15 = NOTcylinder
    24 = NOTsphereUp
    25 = NOTcylinderUp

    PSFmodel: (default = 0)
    0 = GL,
    1 = 3DG,
    2 = hybrid
    3 = GL special
    4 = GL special (for spherical and cylindrical geometry)

    ORDER:         2 = squared, 1 = linear       (default = 1)

    RETURN
    volume of PSF as a function of zexp
    """
    if geo == None:
        geo = 0
    if psfmodel == None:
        psfmodel = 3
    if order == None:
        order = 1

    if order == 1:
        if geo < 10:
            if psfmodel == 0:
                pass
            elif psfmodel == 1:
                pass
            elif psfmodel == 2:
                pass
            elif psfmodel == 3:
                pass
            elif psfmodel == 4:
                pass
        elif 10 <= geo < 15:
            pass
        elif 20 <= geo <= 25:
            pass
        else:
            raise ValueError("geo model is not implemented")

    elif order == 2:
        pass
    else:
        raise ValueError("The value of order is not allowed")

    return (vol=0)







def zscanMultiLayer(z, zoff, parasPSF, model=[], psfmodel="mGL"):
    if model == []:
        raise ValueError("No model is applied.")

    for x in model:

        if x.
