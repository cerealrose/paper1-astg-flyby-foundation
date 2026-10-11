"""Section 4.4 of the paper: osculating hyperbola at perigee for the twelve Earth flybys from JPL Horizons
geocentric state vectors (monopole GM only, as in Anderson et al. 2008), compared with data/earth_flybys.csv.
For the six encounters of Anderson et al. (2008) this checks Vinf and the declinations; for the six later
encounters it is the source of the tabulated geometry.  Needs network access to the Horizons API (runs in
Google Colab); everything else in this repository runs offline.

Frame: geocentric (500@399), ICRF (REF_PLANE='FRAME'), TDB.  The perigee is found on a 1-minute grid and
refined on a 1-second grid."""
from datetime import datetime, timedelta
import numpy as np
import requests
import astg

API = 'https://ssd.jpl.nasa.gov/api/horizons.api'
GM_KM = astg.GM / 1e9
IDS = {'Galileo I': '-77', 'Galileo II': '-77', 'NEAR': '-93', 'Cassini': '-82', 'Rosetta I': '-226', 'MESSENGER': '-236',
       'Rosetta II': '-226', 'EPOXI I': '-140', 'EPOXI II': '-140', 'Rosetta III': '-226', 'EPOXI III': '-140', 'Juno': '-61'}


def parse(result):
    """Rows between $$SOE and $$EOE of a CSV vector table: JD, date, position (km), velocity (km/s)."""
    if '$$SOE' not in result:
        raise RuntimeError(result.strip().splitlines()[-1] if result.strip() else 'empty reply')
    body = result.split('$$SOE')[1].split('$$EOE')[0]
    rows = [l.split(',') for l in body.strip().splitlines() if l.strip()]
    jd = np.array([float(r[0]) for r in rows]); date = [r[1].strip() for r in rows]
    rv = np.array([[float(x) for x in r[2:8]] for r in rows])
    return jd, date, rv[:, :3], rv[:, 3:]


def vectors(sid, start, stop, step):
    p = {'format': 'json', 'COMMAND': f"'{sid}'", 'OBJ_DATA': "'NO'", 'MAKE_EPHEM': "'YES'", 'EPHEM_TYPE': "'VECTORS'",
         'CENTER': "'500@399'", 'REF_PLANE': "'FRAME'", 'REF_SYSTEM': "'ICRF'", 'VEC_TABLE': "'2'", 'OUT_UNITS': "'KM-S'",
         'CSV_FORMAT': "'YES'", 'VEC_LABELS': "'NO'", 'TIME_DIGITS': "'FRACSEC'",
         'START_TIME': f"'{start}'", 'STOP_TIME': f"'{stop}'", 'STEP_SIZE': f"'{step}'"}
    reply = requests.get(API, params=p, timeout=120); reply.raise_for_status()
    return parse(reply.json()['result'])


def osculating(r, v):
    """Vinf, eccentricity, perigee radius, inclination and asymptote declinations (deg) of the osculating hyperbola."""
    h = np.cross(r, v); evec = np.cross(v, h) / GM_KM - r / np.linalg.norm(r)
    e = np.linalg.norm(evec); energy = 0.5 * v @ v - GM_KM / np.linalg.norm(r)
    p_hat = evec / e; h_hat = h / np.linalg.norm(h); q_hat = np.cross(h_hat, p_hat)
    s_in = p_hat / e + np.sqrt(e * e - 1) / e * q_hat          # incoming asymptotic velocity direction
    s_out = -p_hat / e + np.sqrt(e * e - 1) / e * q_hat        # outgoing
    rp = (h @ h / GM_KM) / (1 + e)
    return dict(vinf=np.sqrt(2 * energy), e=e, rp=rp, incl=np.degrees(np.arccos(h_hat[2])),
                din=np.degrees(np.arcsin(s_in[2])), dout=np.degrees(np.arcsin(s_out[2])))


def perigee_state(sid, day):
    d0 = datetime.strptime(day, '%Y-%m-%d'); fmt = lambda t: t.strftime('%Y-%m-%d %H:%M:%S')
    jd, date, r, v = vectors(sid, fmt(d0 - timedelta(days=1)), fmt(d0 + timedelta(days=2)), '1m')
    i = int(np.argmin(np.linalg.norm(r, axis=1)))
    tc = d0 - timedelta(days=1) + timedelta(minutes=i)
    jd, date, r, v = vectors(sid, fmt(tc - timedelta(minutes=3)), fmt(tc + timedelta(minutes=3)), '360')
    j = int(np.argmin(np.linalg.norm(r, axis=1)))
    return date[j], r[j], v[j]


if __name__ == '__main__':
    fl = astg.load_flybys()
    print(f"{'flyby':12s} {'perigee (TDB)':26s} {'hp':>9s} {'Vinf':>7s} {'dVinf(m/s)':>10s} {'d_in':>7s} {'d_out':>7s} {'dcos':>8s} {'tab dcos':>8s} {'incl':>7s}")
    for name, sid in IDS.items():
        rec = fl[name]
        try:
            date, r, v = perigee_state(sid, rec['date'])
        except Exception as ex:
            print(f"{name:12s} query failed: {ex}"); continue
        o = osculating(r, v)
        dc = np.cos(np.radians(o['din'])) - np.cos(np.radians(o['dout']))
        print(f"{name:12s} {date:26s} {o['rp'] - astg.R_EARTH / 1e3:9.1f} {o['vinf']:7.3f} {(o['vinf'] - rec['Vinf_km_s']) * 1e3:+10.1f} "
              f"{o['din']:+7.2f} {o['dout']:+7.2f} {dc:+8.4f} {astg.derived(rec)['dcos']:+8.4f} {o['incl']:7.2f}")
