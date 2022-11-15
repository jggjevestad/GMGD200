# Footpoint latitude
def footlat(a, b, x, lat0):
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = Marc(a, b, lat0) + x

    latf = B/b0 + (3/4*f + 3/8*f**2 + 21/256*f**3)*sin(2*B/b0) \
           + (21/64*f**2 + 21/64*f**3)*sin(4*B/b0) \
           + 151/768*f**3*sin(6*B/b0)
    return latf
