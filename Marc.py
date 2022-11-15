def Marc(a, b, lat):
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = b0*(lat - (3/4*f + 3/8*f**2 + 15/128*f**3)*sin(2*lat)
            + (15/64*f**2 + 15/64*f**3)*sin(4*lat)
            - 35/384*f**3*sin(6*lat))
    return B
