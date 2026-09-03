# Geodetiske funksjoner GMGD200

Oversikt over Python-funksjonshoder basert på formelsamlingen for **GMGD200**.


---

## 1. Vinkelenheter

```python
def deg2rad(deg: float) -> float:
    """Convert degrees to radians."""

def rad2deg(rad: float) -> float:
    """Convert radians to degrees."""

def grad2rad(grad: float) -> float:
    """Convert gradians (gon) to radians."""

def rad2grad(rad: float) -> float:
    """Convert radians to gradians (gon)."""

def dms2deg(d: float, m: float, s: float) -> float:
    """Convert degrees, minutes, seconds to decimal degrees."""

def deg2dms(deg: float) -> tuple[int, int, float]:
    """Convert decimal degrees to (degrees, minutes, seconds)."""

def dms2rad(d: float, m: float, s: float) -> float:
    """Convert degrees, minutes, seconds to radians."""

def rad2dms(rad: float) -> tuple[int, int, float]:
    """Convert radians to (degrees, minutes, seconds)."""
```


---

## 2. Ellipsoide- og krumningsradier

```python
def Mrad(a: float, b: float, lat: float) -> float:
    """Calculate meridional radius of curvature (M)."""

def Nrad(a: float, b: float, lat: float) -> float:
    """Calculate normal radius of curvature (N)."""

def Rm(a: float, b: float, lat: float) -> float:
    """Calculate mean radius of curvature: Rm = sqrt(M * N)."""

def Ra(a: float, b: float, lat: float, az: float) -> float:
    """Calculate radius of curvature for given azimuth (Euler's equation)."""

def Marc(a: float, b: float, lat: float) -> float:
    """Calculate meridional arc distance from equator to latitude."""

def footlat(a: float, b: float, x: float, lat0: float) -> float:
    """Calculate footpoint latitude from arc distance / northing."""
```


---

## 3. Koordinatkonverteringer, ECEF og toposentriske rammer

```python
import numpy as np
from numpy.typing import NDArray

def geod2ECEF(a: float, b: float, lat: float, lon: float, h: float) -> tuple[float, float, float]:
    """Convert geodetic coordinates (lat, lon, h) to ECEF coordinates (X, Y, Z)."""

def ECEF2geod(a: float, b: float, X: float, Y: float, Z: float) -> tuple[float, float, float]:
    """Convert ECEF coordinates (X, Y, Z) to geodetic coordinates (lat, lon, h) using iteration."""

def ECEF2geodb(a: float, b: float, X: float, Y: float, Z: float) -> tuple[float, float, float]:
    """Convert ECEF coordinates (X, Y, Z) to geodetic coordinates using Bowring's closed-form algorithm (1976)."""

def ECEF2geodv(a: float, b: float, X: float, Y: float, Z: float) -> tuple[float, float, float]:
    """Convert ECEF coordinates (X, Y, Z) to geodetic coordinates using Vermeille's method (2004)."""

def ECEF2ned(lat0: float, lon0: float, dX: float, dY: float, dZ: float) -> NDArray:
    """Convert baseline vector (dX, dY, dZ) from ECEF coordinates to NED coordinates."""

def ned2ECEF(lat0: float, lon0: float, n: float, e: float, d: float) -> NDArray:
    """Convert local vector (n, e, d) from NED coordinates to ECEF coordinates."""

def ECEF2enu(lat0: float, lon0: float, dX: float, dY: float, dZ: float) -> NDArray:
    """Convert baseline vector (dX, dY, dZ) from ECEF coordinates to ENU coordinates."""

def enu2ECEF(lat0: float, lon0: float, e: float, n: float, u: float) -> NDArray:
    """Convert local vector (e, n, u) from ENU coordinates to ECEF coordinates."""
```

---

## 4. Geodetiske hovedoppgaver

```python
def geod1(a: float, b: float, lat1: float, lon1: float, az1: float, d: float) -> tuple[float, float, float]:
    """Solve the geodetic direct problem: (lat1, lon1, az1, d) -> (lat2, lon2, az2)."""

def geod2(a: float, b: float, lat1: float, lon1: float, lat2: float, lon2: float) -> tuple[float, float, float]:
    """Solve the geodetic indirect problem: (lat1, lon1, lat2, lon2) -> (lat2_azimuth, az1, d)."""
```

---

## 5. Transversal Mercator og kartprojeksjoner

```python
def geod2TMgrid(
    a: float, b: float, lat: float, lon: float, lat0: float, lon0: float, scale: float, fnorth: float, feast: float
) -> tuple[float, float]:
    """Convert geodetic coordinates (lat, lon) to Transversal Mercator coordinates (north, east)."""

def TMgrid2geod(
    a: float, b: float, north: float, east: float, lat0: float, lon0: float, scale: float, fnorth: float, feast: float
) -> tuple[float, float]:
    """Convert Transversal Mercator coordinates (north, east) to geodetic coordinates (lat, lon)."""

def TMcorr(a: float, b: float, x1: float, y1: float, x2: float, y2: float, lat0: float) -> tuple[float, float]:
    """Calculate Transversal Mercator distance and azimuth corrections (daz, ds)."""

def TMconv1(a: float, b: float, lat: float, lon: float, lon0: float) -> float:
    """Calculate Transversal Mercator meridian convergence gamma from geodetic coordinates (lat, lon)."""

def TMconv2(a: float, b: float, x: float, y: float, lat0: float) -> float:
    """Calculate Transversal Mercator meridian convergence gamma from grid coordinates (x, y)."""

def TMscale1(a: float, b: float, lat: float, lon: float, lon0: float) -> float:
    """Calculate Transversal Mercator point scale factor from geodetic coordinates (lat, lon)."""

def TMscale2(a: float, b: float, x: float, y: float, lat0: float) -> float:
    """Calculate Transversal Mercator point scale factor from grid coordinates (x, y)."""
```

---

## 6. Rotasjoner

```python
def Rx(rx: float) -> NDArray:
    """3D Direction Cosine Matrix for rotation around the x-axis."""

def Ry(ry: float) -> NDArray:
    """3D Direction Cosine Matrix for rotation around the y-axis."""

def Rz(rz: float) -> NDArray:
    """3D Direction Cosine Matrix for rotation around the z-axis."""

def Ce_g(lat: float, lon: float) -> NDArray:
    """Rotate from Earth-Centered Earth-Fixed frame (e-frame) to local geographic frame (g-frame / NED)."""
```

---

## 7. Kovarians og usikkerhet

```python
def std2cov(std_corr: tuple[float, float, float, float, float, float]) -> NDArray:
    """Convert standard deviation and correlation (sigma_x, sigma_y, sigma_z, rho_xy, rho_xz, rho_yz) to 3x3 covariance matrix."""

def cov2std(C: NDArray) -> tuple[float, float, float, float, float, float]:
    """Convert 3x3 covariance matrix to standard deviations and correlation coefficients."""
```

---

## 8. Geometriske reduksjoner og tyngde

```python
def red_dist_chord(d1: float, h1: float, h2: float, R: float, k: float) -> float:
    """Reduce measured chord distance to ellipsoidal arc d4."""

def red_dist_spatial(d2: float, h1: float, h2: float, R: float) -> float:
    """Reduce spatial slope distance d2 to ellipsoidal arc distance d4."""

def red_dist_zenith(d2: float, z: float, h1: float, R: float) -> float:
    """Reduce spatial slope distance d2 to ellipsoidal arc d4 using measured zenith angle z."""

def red_az_skew(lat1: float, lat2: float, az: float, h2: float, a: float, b: float) -> float:
    """Azimuth reduction for skew normals: delta_s."""

def red_az_defl(az: float, z: float, xi: float, eta: float) -> float:
    """Azimuth reduction for deflection of the vertical."""

def zenith_refr_free(d4: float, dH: float, Hm: float, Ra: float) -> float:
    """Calculate refraction-free geometric zenith angle z."""

def grav_height(h: float, M: float, R: float, G: float) -> float:
    """Gravity acceleration as a function of geometric height"""

def norm_grav1980(lat: float) -> float:
    """Normal gravity from the International Gravity Formula 1980 (Somigliana)."""
```

---

