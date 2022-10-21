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

