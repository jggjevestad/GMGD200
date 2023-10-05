# Import
from numpy import pi, arctan2, fmod, fix, sin, cos, tan, arcsin, arctan, sqrt, array


# Modified arctanc (returns quadrant independent angle, e.g. azimuth)
def arctanc(y, x):
    z = arctan2(y, x)
    return fmod(2*pi + z, 2*pi)


# Convert from degree to radian
def deg2rad(deg):
    return deg*(pi/180)


# Convert from radian to degree
def rad2deg(rad):
    return rad*(180/pi)


# Convert from gradian to radian
def grad2rad(grad):
    return grad*(pi/200)


# Convert from radian to gradian
def rad2grad(rad):
    return rad*(200/pi)


# Convert from semicircle to radian
def sc2rad(sc):
    return sc*pi


# Convert from radian to semicircle
def rad2sc(rad):
    return rad/pi


# Convert from degree, minutes, seconds to degree
def dms2deg(dms):
    d = dms[0]
    m = dms[1]
    s = dms[2]

    deg = abs(d) + m/60 + s/3600
    return deg


# Convert from degree to degree, minutes, seconds
def deg2dms(deg):
    frac = abs(deg - int(deg))
    d = int(fix(deg))
    dmin = frac*60

    frac = abs(dmin - int(dmin))
    m = int(fix(dmin))
    s = frac*60
    return [d, m, s]


# Convert from degree, minutes, seconds to radian
def dms2rad(dms):
    deg = dms2deg(dms)
    rad = deg2rad(deg)
    return rad


# Convert from radian to degree, minutes, seconds
def rad2dms(rad):
    deg = rad2deg(rad)
    dms = deg2dms(deg)
    return dms


# Convert from geodic to ECEF on sphere
def geod2ECEFs(R, lat, lon, h):
    
    X = (R + h)*cos(lat)*cos(lon)
    Y = (R + h)*cos(lat)*sin(lon)
    Z = (R + h)*sin(lat)
    
    return X, Y, Z

# Convert from ECEF to geodic on sphere
def ECEF2geods(R, P):
    
    X = P[0]
    Y = P[1]
    Z = P[2]
    
    p = sqrt(X**2 + Y**2)
    lat = arctan(Z/p)
    lon = arctan(Y/X)
    h = p*cos(lat) + Z*sin(lat) - R
    
    return lat, lon, h


# Meridional Radius of Curvature
def Mrad(a, b, lat):
    
    e2 = (a**2 - b**2)/a**2
    M = a*(1 - e2)/(1 - e2*sin(lat)**2)**(3/2)
    
    return M


# Normal Radius of Curvature
def Nrad(a, b, lat):
    
    e2 = (a**2 - b**2)/a**2
    N = a/(1 - e2*sin(lat)**2)**(1/2)
    
    return N


# Convert geodetic coordinates to ECEF coordinates
def geod2ECEF(a, b, lat, lon, h):
        
    N = Nrad(a, b, lat)

    X = (N + h)*cos(lat)*cos(lon)
    Y = (N + h)*cos(lat)*sin(lon)
    Z = ((b**2/a**2)*N + h)*sin(lat)
    
    return X, Y, Z


# Convert ECEF coordinates to geodetic coordinates (iteration)
def ECEF2geod(a, b, X, Y, Z):
   
    e2 = (a**2 - b**2)/a**2

    p = sqrt(X**2 + Y**2)
    lat_new = arctan(Z/p)

    epsilon = 1e-10
    lat = 0

    # for i in range(0, 5):
    while abs(lat_new - lat) > epsilon:
        lat = lat_new
        N = Nrad(a, b, lat)
        lat_new = arctan(Z/p + N*e2*sin(lat)/p)

    lat = lat_new
    lon = arctan(Y/X)
    h = p*cos(lat) + Z*sin(lat) - N*(1 - e2*sin(lat)**2)

    return lat, lon, h


# Convert ECEF coordinates to geodetic coordinates (Bowring)
def ECEF2geodb(a, b, X, Y, Z):
    
    e2 = (a**2 - b**2)/a**2
    e2m = (a**2 - b**2)/b**2

    p = sqrt(X**2 + Y**2)
    mu = arctan(Z*a/(p*b))

    lat = arctan((Z + e2m*b*sin(mu)**3)/(p - e2*a*cos(mu)**3))
    lon = arctan(Y/X)
    h = p*cos(lat) + Z*sin(lat) - Nrad(a, b, lat)*(1 - e2*sin(lat)**2)

    return lat, lon, h


# Convert ECEF coordinates to ENU coordinates
def ECEF2enu(lat0, lon0, dX, dY, dZ):
    
    P1 = array([[dX],
                [dY],
                [dZ]])
    
    Ce_g = array([[-sin(lon0), cos(lon0), 0],
                  [-sin(lat0)*cos(lon0), -sin(lat0)*sin(lon0), cos(lat0)],
                  [cos(lat0)*cos(lon0), cos(lat0)*sin(lon0), sin(lat0)]])
    
    P2 = Ce_g@P1
    
    return P2[0, 0], P2[1, 0], P2[2, 0]


def Rx(alpha):
    return array([[1 ,0 , 0],
                  [0, cos(alpha), sin(alpha)],
                  [0, -sin(alpha), cos(alpha)]])

def Ry(beta):
    return array([[cos(beta) ,0 , -sin(beta)],
                  [0, 1, 0],
                  [sin(beta), 0, cos(beta)]])

def Rz(gamma):
    return array([[cos(gamma) ,sin(gamma) , 0],
                  [-sin(gamma), cos(gamma), 0],
                  [0, 0, 1]])


# Convert geodetic coordinates to n-vector
def geod2n(lat, lon):
    return array([[cos(lat)*cos(lon)],
                  [cos(lat)*sin(lon)],
                  [sin(lat)]])


# Convert n-vector to geodetic coordinates
def n2geod(n):
    
    lat = arctan2(n[2, 0], sqrt(n[0, 0]**2 + n[1, 0]**2))
    lon = arctan2(n[1, 0], n[0, 0])
    
    return lat, lon


# Geodetic direct problem
def geod1(a, b, lat1, lon1, az1, d):
    
    f = (a - b)/a
    e2m = (a**2 - b**2)/b**2

    beta1 = arctan(b/a*tan(lat1))
    az0 = arcsin(sin(az1)*cos(beta1))
    sigma1 = arctan(tan(beta1)/cos(az1))

    g = e2m*cos(az0)**2
    H = 1/8*g - 1/16*g**2 + 37/1024*g**3
    b0 = b*(1 + 1/4*g - 3/64*g**2 + 5/256*g**3)

    d1 = b0*(sigma1 - H*sin(2*sigma1) - H**2/4*sin(4*sigma1) - H**3/6*sin(6*sigma1))
    d2 = d1 + d

    sigma2 = d2/b0 + (H - 3/4*H**3)*sin(2*d2/b0) + 5/4*H**2*sin(
        4*d2/b0) + 29/12*H**3*sin(6*d2/b0)
    sigma = sigma2 - sigma1

    X = cos(beta1)*cos(sigma) - sin(beta1)*sin(sigma)*cos(az1)
    Y = sin(sigma)*sin(az1)
    Z = sin(beta1)*cos(sigma) + cos(beta1)*sin(sigma)*cos(az1)

    beta2 = arctan(Z/sqrt(X**2 + Y**2))
    dlon = arctan(Y/X)

    K = (f + f**2)/4*cos(az0)**2 - f**2/4*cos(az0)**4
    dlon = dlon - f*sin(az0)*((1 - K - K**2)*sigma + K*sin(sigma)*cos(sigma1 + sigma2)
                              + K**2*sin(sigma)*cos(sigma)*cos(2*(sigma1 + sigma2)))

    lat2 = arctan(a/b*tan(beta2))
    lon2 = lon1 + dlon
    az2 = arctanc(sin(az1)*cos(beta1), (cos(beta1)*cos(sigma)*cos(az1) - 
sin(beta1)*sin(sigma)))

    if az2 < pi:
        az2 = az2 + pi
    else:
        az2 = az2 - pi

    return lat2, lon2, az2


# Geodetic indirect problem
def geod2(a, b, lat1, lon1, lat2, lon2):

    f = (a - b)/a
    e2m = (a**2 - b**2)/b**2

    beta1 = arctan(b/a*tan(lat1))
    beta2 = arctan(b/a*tan(lat2))

    epsilon = 1e-10
    dlon_new = lon2 - lon1
    dlon = 0

    while abs(dlon_new - dlon) > epsilon:
        dlon = dlon_new

        X = cos(beta1)*sin(beta2) - sin(beta1)*cos(beta2)*cos(dlon)
        Y = cos(beta2)*sin(dlon)
        Z = sin(beta1)*sin(beta2) + cos(beta1)*cos(beta2)*cos(dlon)

        sigma = arctan(sqrt(X**2 + Y**2)/Z)
        az1 = arctanc(Y, X)
        az0 = arcsin(sin(az1)*cos(beta1))

        sigma1 = arctan(tan(beta1)/cos(az1))
        sigma2 = sigma1 + sigma

        K = (f + f**2)/4*cos(az0)**2 - f**2/4*cos(az0)**4

        dlon_new = (lon2 - lon1) + f*sin(az0)*((1 - K - K**2)*sigma + K*sin(sigma)*cos(sigma1 + sigma2) 
                                               + K**2*sin(sigma)*cos(sigma)*cos(2*(sigma1 + sigma2)))

    dlon = dlon_new
    az2 = arctanc(cos(beta1)*sin(dlon), (cos(beta1)*sin(beta2)*cos(dlon) - 
sin(beta1)*cos(beta2)))

    if az2 < pi:
        az2 = az2 + pi
    else:
        az2 = az2 - pi

    g = e2m*cos(az0)**2
    H = 1/8*g - 1/16*g**2 + 37/1024*g**3
    b0 = b*(1 + 1/4*g - 3/64*g**2 + 5/256*g**3)

    d = b0*(sigma - 2*H*sin(sigma)*cos(sigma1 + sigma2)
            - H**2/2*sin(2*sigma)*cos(2*(sigma1 + sigma2))
            - H**3/3*sin(3*sigma)*cos(3*(sigma1 + sigma2)))

    return az1, az2, d


# Footpoint latitude
def footlat(a, b, x, lat0):
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = Marc(a, b, lat0) + x

    latf = B/b0 + (3/4*f + 3/8*f**2 + 21/256*f**3)*sin(2*B/b0) \
           + (21/64*f**2 + 21/64*f**3)*sin(4*B/b0) \
           + 151/768*f**3*sin(6*B/b0)
    return latf


# Meridional arc distance
def Marc(a, b, lat):
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = b0*(lat - (3/4*f + 3/8*f**2 + 15/128*f**3)*sin(2*lat)
            + (15/64*f**2 + 15/64*f**3)*sin(4*lat)
            - 35/384*f**3*sin(6*lat))
    return B


# Convert geodetic coordinates to Transversal Mercator coordinates
def geod2TMgrid(a, b, lat, lon, lat0, lon0, scale, fnorth, feast):
    B = Marc(a, b, lat) - Marc(a, b, lat0)
    N = Nrad(a, b, lat)
    e2 = (a**2 - b**2)/a**2
    eps2 = e2/(1 - e2)*cos(lat)**2
    dlon = lon - lon0

    x = B + 1/2*dlon**2*N*sin(lat)*cos(lat) \
        + 1/24*dlon**4*N*sin(lat)*cos(lat)**3*(5 - tan(lat)**2 + 9*eps2 + 4*eps2**2) \
        + 1/720*dlon**6*N*sin(lat)*cos(lat)**5*(61 - 58*tan(lat)**2 + tan(lat)**4)

    y = dlon*N*cos(lat) + 1/6*dlon**3*N*cos(lat)**3*(1 - tan(lat)**2 + eps2) \
        + 1/120*dlon**5*N*cos(lat)**5*(5 - 18*tan(lat)**2 + tan(lat)**4)

    north = x*scale
    east = y*scale

    north = north + fnorth
    east = east + feast

    return north, east


# Convert Transversal Mercator coordinates to geodetic coordinates
def TMgrid2geod(a, b, north, east, lat0, lon0, scale, fnorth, feast):
    north = north - fnorth
    east = east - feast

    x = north/scale
    y = east/scale

    latf = footlat(a, b, x, lat0)

    e2 = (a**2 - b**2)/a**2

    Mf = Mrad(a, b, latf)
    Nf = Nrad(a, b, latf)
    eps2f = e2/(1 - e2)*cos(latf)**2

    lat = latf - 1/2*y**2*tan(latf)/(Mf*Nf) \
          + 1/24*y**4*tan(latf)/(Mf*Nf**3)*(5 + 3*tan(latf)**2 + eps2f - 9*eps2f*tan(latf)**2 
- 4*eps2f**2) \
          - 1/720*y**6*tan(latf)/(Mf*Nf**5)*(61 + 90*tan(latf)**2 + 45*tan(latf)**4)

    dlon = y/(Nf*cos(latf)) \
           - 1/6*y**3/(Nf**3*cos(latf))*(1 + 2*tan(latf)**2 + eps2f) \
           + 1/120*y**5/(Nf**5*cos(latf))*(5 + 28*tan(latf)**2 + 24*tan(latf)**4)

    lon = dlon + lon0

    return lat, lon
