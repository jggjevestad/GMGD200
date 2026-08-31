# GMGD200 – Geodesi (NMBU)

## 1. Course & Project Overview
**GMGD200 (Geodesi)** is an academic course in Geodesy offered at the Norwegian University of Life Sciences (NMBU). The course covers geometric and physical geodesy, reference frames, map projections, coordinate transformations, spherical trigonometry, GNSS vector reductions, and error propagation.

---

## 2. Core Reference Frames & Datums

### Horizontal & 3D Datums
* **EUREF89**:
  * Official 3D/horizontal reference frame in Norway (static European realization tied to ITRF/ETRF89).
  * Reference Ellipsoid: **GRS80** ($a = 6\,378\,137.0\text{ m}$, $f = 1/298.257222101$).
* **NGO1948**:
  * Historic Norwegian datum based on the Oslo Observatory Datum (1844).
  * Reference Ellipsoid: **Bessel 1841** ($a = 6\,377\,397.155\text{ m}$, $1/f = 299.1528128$).
  * Grid: 8 Gauss-Krüger projection zones (NGO zones I–VIII).
* **ED50 (European Datum 1950)**:
  * Post-war international datum (Hayford / International 1924 ellipsoid: $a = 6\,378\,388.0\text{ m}$, $1/f = 297.0$).
* **WGS84 / ITRF**:
  * Global geocentric dynamic systems used by GNSS / satellite positioning.

### Vertical Reference Frames
* **NN2000**:
  * Official national vertical reference frame in Norway (realization of EVRF2007).
  * Uses **Normal heights ($H^N$)** relative to the quasigeoid.
* **NN1954**:
  * Historic vertical datum based on mean sea level in Tregde (orthometric heights).
* **Fundamental relationship**:
  $$h = H + N$$
  *(where $h$ is ellipsoidal height, $H$ is physical/normal height, and $N$ is geoid/quasigeoid height).*

---

## 3. Map Projections

### UTM (Universal Transverse Mercator)
* Conformal Gauss-Krüger projection with 6° zone width.
* Standard zones in Norway: **UTM 32 (9°E)**, **UTM 33 (15°E)**, and **UTM 35 (27°E)**.
* Scale factor at central meridian: $k_0 = 0.9996$ (max $-400\text{ ppm}$ distortion).
* False Easting: $500\,000\text{ m}$, False Northing: $0\text{ m}$.

### NTM (Norsk Transversal Mercator)
* Secondary official projection designed for the construction and engineering industry.
* 1° zone widths (zones NTM 5 to NTM 30, where zone number equals central meridian in °E).
* Scale factor at central meridian: $k_0 = 1.0$ (virtually zero distortion along central meridian).
* False Easting: $100\,000\text{ m}$ (or $1\,000\,000 + \text{zone}\cdot 100\,000\text{ m}$).

---

## 4. Key Mathematical & Geodetic Concepts

1. **Coordinate Conversions**:
   * Cartesian geocentric $(X, Y, Z) \leftrightarrow$ Ellipsoidal geodetic $(\phi, \lambda, h)$.
   * Conversion methods: Iterative method vs. **Bowring’s closed-form algorithm**.
2. **Transformations**:
   * 7-parameter 3D Helmert transformation ($T_X, T_Y, T_Z, D, \alpha, \beta, \gamma$).
   * 2D conformal Helmert transformations.
   * TIN/grid transformations (e.g. SKTrans / Proj grid models for NGO1948 $\leftrightarrow$ EUREF89).
3. **Spherical Trigonometry & Navigation**:
   * Great circle navigation (orthodromes), distance, and forward/backward azimuths.
   * Loxodromes (rhumb lines).
   * Intersecting great circles (e.g., flight path intersections).
4. **Physical Geodesy & Deflection of the Vertical**:
   * Astrogeodetic coordinates $(\Phi, \Lambda)$ vs. Geodetic coordinates $(\phi, \lambda)$.
   * Deflections of the vertical components:
     $$\xi = \Phi - \phi, \quad \eta = (\Lambda - \lambda)\cos\phi$$
5. **GNSS Baseline Reductions**:
   * Converting baseline vectors $(\Delta X, \Delta Y, \Delta Z)$ to topocentric coordinates (azimuth, zenith distance, slope distance).
   * Covariance error propagation and mapping to map grid coordinates.

---

## 5. Exercise Structure (Øving 1–10)

| Exercise | Topic / Focus |
| :--- | :--- |
| **Øving 01** | Ellipsoidal geometry, $(\phi, \lambda, h) \leftrightarrow (X, Y, Z)$, Bowring's method vs. iteration. |
| **Øving 02** | GNSS baseline vectors, topocentric reduction (azimuth, distance, $\Delta h$), datum shifts. |
| **Øving 03** | Spherical Earth approximations, chord height, arc distance calculations. |
| **Øving 04** | Deflection of the vertical, Greenwich 0-meridian physical shift, Airy 1830 vs. WGS84. |
| **Øving 05** | 2D/3D coordinate transformations and parameter estimation. |
| **Øving 06** | Spherical trigonometry, great circle navigation (routes, azimuths, intersections). |
| **Øving 07** | Differential map projection equations, Mercator conformal vs. equal-area projection in Python. |
| **Øving 08** | Ellipsoidal line reduction and transformation to EUREF89. |
| **Øving 09** | Coordinate transformation parameters between NGO1948, ED50, and EUREF89 (Helmert). |
| **Øving 10** | GNSS baseline to UTM map grid projection and covariance error propagation. |

---

## 6. Development & Computational Tools
* **Language**: Python 3
* **Libraries**: `numpy`, `scipy`, `matplotlib`, `pyproj`
* **Reference Specifications**: Kartverket reference frame guides, EPSG registry standards
