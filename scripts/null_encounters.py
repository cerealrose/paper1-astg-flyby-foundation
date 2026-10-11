"""Table 8 of the paper: the seven null encounters against the ASTG law (Eq. 12, lambda1 from Eq. 16) and the
Anderson relation converted to perigee, dVp = kappa_A dcos Vinf^2 / Vp.  Deviations are in units of the combined
uncertainty: for the law, sqrt(sigma_obs^2 + (prediction * se/lambda1)^2); for the Anderson relation, sigma_obs."""
import numpy as np
import astg


def table8(flybys=None):
    flybys = flybys or astg.load_flybys()
    fit = astg.anomalous_fit(flybys)
    lam, se = fit['lam'], fit['se']
    rows = []
    for k in astg.NULLS:
        f = flybys[k]
        a = astg.astg_dvp(f, lam); s_a = abs(a) * se / lam
        an = astg.anderson_dvp(f)
        o, s = f['dVp_mm_s'], f['sigma_dVp_mm_s']
        rows.append(dict(spacecraft=k, observed=o, sigma=s, astg=a, anderson=an,
                         z_astg=abs(a - o) / np.hypot(s, s_a), z_anderson=abs(an - o) / s))
    return rows, lam, se


if __name__ == '__main__':
    rows, lam, se = table8()
    print('Table 8. Null encounters (mm/s; deviations in combined standard deviations)')
    print(f"{'Spacecraft':12s} {'observed':>16s} {'ASTG':>8s} {'Anderson':>9s} {'ASTG (s)':>9s} {'Anderson (s)':>13s}")
    for r in rows:
        print(f"{r['spacecraft']:12s} {r['observed']:+8.4f} +- {r['sigma']:.3f} {r['astg']:+8.3f} {r['anderson']:+9.3f} {r['z_astg']:9.1f} {r['z_anderson']:13.1f}")
    print(f"chi^2 over the seven nulls: ASTG {sum(r['z_astg']**2 for r in rows):.0f}, Anderson {sum(r['z_anderson']**2 for r in rows):.0f}")
    better = sum(r['z_astg'] < r['z_anderson'] for r in rows)
    print(f"the ASTG law is closer than the Anderson relation for {better} of the seven")
    juno = [r for r in rows if r['spacecraft'] == 'Juno'][0]
    s_det = 0.0057      # 1-sigma detection limit of the Juno tracking (companion paper, Section 4.1)
    print(f"Juno at the detection limit {s_det} mm/s: {juno['astg']/np.hypot(s_det, juno['astg']*se/lam):.1f} sigma")
