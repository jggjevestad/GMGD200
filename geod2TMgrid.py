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
