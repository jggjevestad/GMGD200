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

# Footpoint latitude
def footlat(a, b, x, lat0):
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = Marc(a, b, lat0) + x

    latf = B/b0 + (3/4*f + 3/8*f**2 + 21/256*f**3)*sin(2*B/b0) \
           + (21/64*f**2 + 21/64*f**3)*sin(4*B/b0) \
           + 151/768*f**3*sin(6*B/b0)
    return latf
