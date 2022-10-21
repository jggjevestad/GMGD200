# Geodetic indirect problem
def geod2(a, b, lat1, lon1, lat2, lon2):
    az0 = None
    az1 = None
    sigma = None
    sigma1 = None
    sigma2 = None

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

        dlon_new = dlon + f*sin(az0)*((1 - K - K**2)*sigma + K*sin(sigma)*cos(sigma1 + sigma2)
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

