# Modified arctanc (returns quadrant independent angle, e.g. azimuth)
def arctanc(y, x):
    z = arctan2(y, x)
    return fmod(2*pi + z, 2*pi)
