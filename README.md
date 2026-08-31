# GMGD200 – Geodesi (NMBU)

Dette repositoriet inneholder øvinger, beregningsskript og kursmateriell for emnet **GMGD200 (Geodesi)** ved **Norges miljø- og biovitenskapelige universitet (NMBU)**.

---

## Emneoversikt

**GMGD200 (Geodesi)** dekker grunnleggende og videregående prinsipper innen geometrisk og fysisk geodesi, referanserammer, koordinatsystemer, kartprojeksjoner og reduksjon av satellittbaserte posisjonsmålinger (GNSS). Emnet vektlegger både det teoretiske fundamentet og praktiske numeriske beregninger ved hjelp av Python.

### Sentrale temaer
- **Geometrisk og fysisk geodesi**: Jordellipsoider, geoide- og kvasigeoidemodeller, astronomiske koordinater og loddavvik ($\xi, \eta$).
- **Referanserammer og datum**:
  - **Horisontale / 3D-systemer**: EUREF89 (GRS80-ellipsoiden), NGO1948 (Bessel 1841-ellipsoiden), ED50 (Hayford / Int. 1924) samt globale WGS84 og ITRF.
  - **Høydereferansesystemer / Vertikale datum**: NN2000 (normalhøyder) og historiske NN1954 (ortometriske høyder).
  - Grunnleggende høyderelasjon: $h = H + N$.
- **Kartprojeksjoner**:
  - **UTM (Universal Transverse Mercator)**: Sonene 32N, 33N og 35N i Norge ($k_0 = 0.9996$).
  - **NTM (Norsk Transversal Mercator)**: $1^\circ$-brede projeksjonssoner (NTM 5–30) tilpasset bygg- og anleggssektoren med minimal målestokkforvrengning langs sentralmeridianen ($k_0 = 1.0$).
- **Koordinattransformasjoner og geodetiske algoritmer**:
  - Lukket form (Bowrings metode) og iterative koordinatkonverteringer: $(\phi, \lambda, h) \leftrightarrow (X, Y, Z)$.
  - 2D- og 3D-Helmert-transformasjoner (7-parameter).
  - Sfærisk trigonometri, storsirkelnavigasjon (ortodromer), loksodromer og skjæringspunkter mellom ruter.
- **GNSS og feilforplantning**:
  - GNSS-vektorreduksjon: Konvertering fra geosentriske komponenter $(\Delta X, \Delta Y, \Delta Z)$ til lokale toposentriske koordinater (asimut, senitavstand, skråavstand).
  - Kovariansmatrisetransformasjoner og feilforplantning til kartplanet.

---

## Øvingsoversikt (Øving 1–10)

| Øving | Tema / Fokus |
| :--- | :--- |
| **Øving 01** | Ellipsoidegeometri, $(\phi, \lambda, h) \leftrightarrow (X, Y, Z)$, Bowrings metode vs. iterativ metode |
| **Øving 02** | GNSS-vektorer, toposentrisk reduksjon (asimut, avstand, $\Delta h$), datumskift |
| **Øving 03** | Sfæriske jordtilnærminger, kordehøyde, bueavstandsberegninger |
| **Øving 04** | Loddavvik, fysisk forskyvning av Greenwich-nullmeridianen, Airy 1830 vs. WGS84 |
| **Øving 05** | 2D- og 3D-koordinattransformasjoner og parameterestimering |
| **Øving 06** | Sfærisk trigonometri, storsirkelnavigasjon (ruter, asimuter, skjæringspunkter) |
| **Øving 07** | Differensielle kartprojeksjonsligninger, konform Mercator vs. arealtro projeksjoner |
| **Øving 08** | Ellipsoidisk linjereduksjon og transformasjon til EUREF89 |
| **Øving 09** | Koordinattransformasjonsparametere mellom NGO1948, ED50 og EUREF89 (Helmert) |
| **Øving 10** | GNSS-vektor til UTM-kartprojeksjon og kovariansfeilforplantning |

---

## Beregningsverktøy

- **Programmeringsspråk**: Python 3
- **Hovedbiblioteker**: `numpy`, `scipy`, `matplotlib`, `pyproj`
- **Geodetiske standarder**: Kartverkets standarder og veiledere for referanserammer samt EPSG-registeret