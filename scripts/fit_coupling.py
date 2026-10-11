"""Determination of the coupling (Section 5 of the paper):
  Table 5  single-encounter couplings, Eq. (13), with sigma propagated from sigma(dVp)
  Eq. (16) no-intercept least-squares fit to the four anomalous encounters, R^2, t, p, and d = lambda1 GM / c0^2
  Table 6  predictions with 95% prediction intervals dVp +- t_{0.025,3} sigma_res Vp sqrt(1 + h)
  Table 7  leverage, standardised residual, Cook's distance and leave-one-out refits; PRESS and RSS
  Eq. (17) control fit with Galileo II added."""
import numpy as np
from scipy import stats
import astg

ANOM = astg.ANOMALOUS


def single_encounter_couplings(flybys):
    rows = []
    for k in ['Galileo I', 'NEAR', 'Cassini', 'Rosetta I', 'Galileo II', 'MESSENGER']:
        f = flybys[k]; d = astg.derived(f)
        lam_i = f['dVp_mm_s'] * 1e-3 / d['vp'] / d['X']
        sig_i = f['sigma_dVp_mm_s'] * 1e-3 / d['vp'] / abs(d['X'])
        rows.append(dict(spacecraft=k, lam_1e3=lam_i / 1e3, sigma_1e3=sig_i / 1e3))
    return rows


def main_fit(flybys):
    fit = astg.anomalous_fit(flybys)
    fit['d_m'] = fit['lam'] * astg.GM / astg.C0**2
    fit['sigma_d_m'] = fit['se'] * astg.GM / astg.C0**2
    fit['d_over_rp_NEAR'] = fit['d_m'] / astg.derived(flybys['NEAR'])['rp']
    return fit


def predictions(flybys, fit):
    x = np.array([astg.derived(flybys[k])['X'] for k in ANOM])
    t975 = stats.t.ppf(0.975, fit['dof']); sres = np.sqrt(fit['s2'])
    rows = []
    for k in ANOM + ['Galileo II']:
        f = flybys[k]; d = astg.derived(f)
        pred = fit['lam'] * d['X'] * d['vp'] * 1e3
        h = d['X']**2 / (x @ x)
        half = t975 * sres * np.sqrt(1 + h) * d['vp'] * 1e3
        rows.append(dict(spacecraft=k, observed=f['dVp_mm_s'], sigma=f['sigma_dVp_mm_s'], astg=pred, lo=pred - half, hi=pred + half))
    return rows, t975


def leave_one_out(flybys, fit):
    x = np.array([astg.derived(flybys[k])['X'] for k in ANOM])
    y = np.array([flybys[k]['dVp_mm_s'] * 1e-3 / astg.derived(flybys[k])['vp'] for k in ANOM])
    vp = np.array([astg.derived(flybys[k])['vp'] for k in ANOM])
    h = fit['leverage']; s2 = fit['s2']; res = fit['res']
    sres = res / np.sqrt(s2 * (1 - h))
    cook = sres**2 * h / (1 - h)            # one fitted parameter
    rows = []
    for i, k in enumerate(ANOM):
        m = np.arange(len(ANOM)) != i
        lam_i = x[m] @ y[m] / (x[m] @ x[m])
        rows.append(dict(encounter=k, leverage=h[i], std_residual=sres[i], cooks_d=cook[i], lam_loo=lam_i,
                         shift_percent=(lam_i / fit['lam'] - 1) * 100))
    press = np.sum(((res / (1 - h)) * vp * 1e3)**2)       # (mm/s)^2
    rss = np.sum((res * vp * 1e3)**2)
    return rows, press, rss


def control_fit(flybys):
    return astg.anomalous_fit(flybys, names=ANOM + ['Galileo II'])


if __name__ == '__main__':
    fl = astg.load_flybys()
    print('Table 5. Single-encounter couplings, Eq. (13)')
    for r in single_encounter_couplings(fl):
        print(f"  {r['spacecraft']:11s} lambda1 = {r['lam_1e3']:7.3f} e3   sigma = {r['sigma_1e3']:7.3f} e3")
    fit = main_fit(fl)
    print(f"\nEq. (16): lambda1 = ({fit['lam']/1e3:.3f} +- {fit['se']/1e3:.3f}) e3   [unrounded {fit['lam']:.1f} +- {fit['se']:.1f}]")
    print(f"  R^2 = {fit['r2_mean']:.3f} about the mean ({fit['r2_origin']:.3f} about the origin), t = {fit['t']:.3f}, p = {fit['p']:.4f}, dof = {fit['dof']}")
    print(f"  axial displacement d = lambda1 GM/c0^2 = {fit['d_m']:.1f} +- {fit['sigma_d_m']:.1f} m;  d/rp (NEAR) = {fit['d_over_rp_NEAR']:.2e}")
    rows, t975 = predictions(fl, fit)
    print(f"\nTable 6. Predictions with 95% prediction intervals (t_0.025,3 = {t975:.2f}), mm/s")
    for r in rows:
        print(f"  {r['spacecraft']:11s} observed {r['observed']:+7.3f} +- {r['sigma']:.3f}   ASTG {r['astg']:+7.3f}   [{r['lo']:+.2f}, {r['hi']:+.2f}]")
    loo, press, rss = leave_one_out(fl, fit)
    print('\nTable 7. Leave-one-out sensitivity')
    for r in loo:
        print(f"  {r['encounter']:10s} h = {r['leverage']:.4f}  std. residual = {r['std_residual']:+.2f}  Cook's D = {r['cooks_d']:.3f}"
              f"  lambda1(-i) = {r['lam_loo']:.0f} ({r['shift_percent']:+.1f}%)")
    print(f"  PRESS = {press:.1f} (mm/s)^2, in-sample RSS = {rss:.1f} (mm/s)^2, ratio {press/rss:.2f}")
    c = control_fit(fl)
    print(f"\nEq. (17), Galileo II added: lambda1 = ({c['lam']/1e3:.3f} +- {c['se']/1e3:.3f}) e3, R^2 = {c['r2_mean']:.3f}, t = {c['t']:.2f}, p = {c['p']:.4f}, n = {c['dof']+1}")
