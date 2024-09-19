# Konverter fra ECEF koordinater til geodetiske koordinater (ellipsoide)
def ECEF2geod(a, b, X, Y, Z, tolerance=1e-10):
    e2 = (a**2 - b**2) / a**2
    rho = sqrt(X**2 + Y**2)
    
    # Initialiserer
    converged = False
    lat_new = arctan(Z / rho)
    lat = 0

    while abs(lat_new - lat) > tolerance:
        lat = lat_new
        N = Nrad(a, b, lat)
        lat_new = arctan(Z / rho + N * e2 * sin(lat) / rho)

    lat = lat_new
    lon = arctan(Y / X)
    h = rho * cos(lat) + Z * sin(lat) - N * (1 - e2 * sin(lat)**2)
    
    return lat, lon, h