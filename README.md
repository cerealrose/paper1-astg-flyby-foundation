# Application of the ASTG-Model to Earth Flyby Encounters: Plausible Origins and Velocity Anomalies

P. A. Nesvinga<sup>1,2</sup>, B. C. Jones<sup>3</sup> and G. G. Nyambuya<sup>1,2</sup>

<sup>1</sup> Fundamental Theoretical Physics Group, National University of Science and Technology, Bulawayo, Zimbabwe<br>
<sup>2</sup> Department of Applied Physics, National University of Science and Technology, Bulawayo, Zimbabwe<br>
<sup>3</sup> Department of Applied Mathematics, National University of Science and Technology, Bulawayo, Zimbabwe

Data and code for the first paper of the series, submitted to *Advances in Space Research* (2026). The paper fits the
flyby law of the Azimuthally Symmetric Theory of Gravitation (ASTG),

ΔV<sub>p</sub>/V<sub>p</sub> = λ<sub>1</sub> R<sub>g</sub> Δcos δ / [2 r<sub>p</sub>(1 + ε)], with R<sub>g</sub> = 2GM<sub>⊕</sub>/c<sub>0</sub><sup>2</sup>,

to the four anomalous Earth flybys, finds λ<sub>1</sub> = (3.902 ± 0.576) × 10<sup>3</sup>, and tests it against the seven null encounters.

## What the paper says is here

Appendix C of the paper reads:

> The flyby dataset of Tables 2 and 3, the derived quantities and fits of Tables 4–8, the exchange budget of Table 1 and the code producing Figures 1 and 2 are available at https://github.com/cerealrose.

| Item in the paper | Where it is |
|---|---|
| Flyby dataset, Tables 2 and 3 | `data/earth_flybys.csv` |
| Exchange budget, Table 1 | `scripts/table1_exchange_budget.py` |
| Derived orbital quantities, Table 4 | `scripts/table4_derived_quantities.py` |
| Single-encounter couplings, Table 5; global fit, Eq. (16); predictions, Table 6; leave-one-out sensitivity, Table 7; control fit, Eq. (17) | `scripts/fit_coupling.py` |
| Null encounters against the law and the Anderson relation, Table 8 | `scripts/null_encounters.py` |
| Figures 1 and 2 | `scripts/figures.py`, written to `figures/` |
| Geometry of the six later encounters from JPL Horizons, Section 4.4 | `scripts/horizons_check.py` (needs network access) |
| Every tabulated number against the printed value | `scripts/check_against_paper.py` |

## Running it

Python 3.9 or later with the packages in `requirements.txt`:

```
pip install -r requirements.txt
cd scripts
python check_against_paper.py
```

`check_against_paper.py` recomputes the numbers of Tables 1 and 4 to 8, Eqs. (16) and (17) and the χ<sup>2</sup> sums of Section 6,
compares each with the value printed in the paper at the printed precision, and prints one PASS or FAIL line per number
(161 in all). Each of the other scripts prints its own table. `python figures.py` redraws both figures. Everything runs
offline in a few seconds except `horizons_check.py`, which queries the JPL Horizons API and also runs in Google Colab.

## Data

`data/earth_flybys.csv` has one row per documented Earth flyby, twelve in all.

| Column | Meaning |
|---|---|
| `spacecraft`, `date` | Encounter and date of closest approach |
| `class` | `anomalous` (the four encounters that fix the coupling), `null`, or `control` (Galileo II) |
| `dVp_mm_s`, `sigma_dVp_mm_s` | Perigee velocity anomaly and its uncertainty, mm/s (Table 2) |
| `anomaly_source` | Source of the anomaly: Jouannic et al. (2015), or Anderson et al. (2008) scaled to perigee by V<sub>∞</sub>/V<sub>p</sub> |
| `dVinf_mm_s`, `sigma_dVinf_mm_s` | Change of hyperbolic excess speed of Anderson et al. (2008), used for Table 1 |
| `hp_km` | Perigee altitude above the 6371 km mean radius (Table 3) |
| `Vinf_km_s` | Hyperbolic excess speed (Table 3) |
| `delta_in_deg`, `delta_out_deg` | Declinations of the incoming and outgoing asymptotic velocity vectors (Table 3) |
| `geometry_source` | Anderson et al. (2008), or the JPL Horizons osculating state at perigee |

Constants (in `scripts/astg.py`): GM<sub>⊕</sub> = 3.986004418 × 10<sup>14</sup> m<sup>3</sup> s<sup>−2</sup> (IERS Conventions 2010), R<sub>⊕</sub> = 6371 km,
ω<sub>⊕</sub> = 7.292115 × 10<sup>−5</sup> rad s<sup>−1</sup>, c<sub>0</sub> = 299 792 458 m s<sup>−1</sup>, so that R<sub>g</sub> = 8.870 × 10<sup>−3</sup> m and the Anderson
constant κ<sub>A</sub> = 2ω<sub>⊕</sub>R<sub>⊕</sub>/c<sub>0</sub> = 3.099 × 10<sup>−6</sup>.

## Companion papers

| Paper | Repository |
|---|---|
| A Drag Model for the Low-Perigee Galileo II Encounter and the Manoeuvre-Affected and Null Earth Gravity Assists | [paper2-astg-flyby-anomalies](https://github.com/cerealrose/paper2-astg-flyby-anomalies) |
| Application of the ASTG-Model to Earth Flyby Encounters: On a Pre-Flyby Predictive Test of the ASTG | [paper3-astg-flyby-prediction](https://github.com/cerealrose/paper3-astg-flyby-prediction) |

The three repositories share `data/earth_flybys.csv` and `scripts/astg.py`.

## Earlier version

The derivation checks that accompanied an earlier version of the manuscript are kept in the history of this repository
under the tag `v1`.

## Licence and citation

Code is under the MIT licence and the dataset under CC BY 4.0; see `LICENSE`. Citation details are in `CITATION.cff`.

Contact: P. A. Nesvinga, praise.nesvinga@nust.ac.zw
