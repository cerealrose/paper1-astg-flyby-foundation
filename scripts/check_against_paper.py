"""Recompute every tabulated number of the paper and compare it with the printed value at the printed precision.
Prints one line per check and a summary; exits with status 1 if any check fails."""
import sys
import astg
from table1_exchange_budget import table1
from table4_derived_quantities import table4
from fit_coupling import single_encounter_couplings, main_fit, predictions, leave_one_out, control_fit
from null_encounters import table8

results = []


def check(label, value, printed, decimals):
    ok = round(value, decimals) == round(printed, decimals) or abs(value - printed) <= 0.5 * 10**(-decimals) + 1e-12
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {label:55s} computed {value:+.{decimals + 3}f}  printed {printed:+.{decimals}f}")


fl = astg.load_flybys()

# Table 1: dE, dLz, sigma, dLz/J
T1 = {'Galileo I': (35.1, 4.81, 0.37, 4.8), 'NEAR': (92.2, 12.65, 0.01, 14.4), 'Cassini': (-32.0, -4.39, 2.20, -3.1),
      'Rosetta I': (6.95, 0.95, 0.02, 1.1), 'Galileo II': (-40.8, -5.60, 1.22, -6.0), 'MESSENGER': (0.08, 0.011, 0.006, 0.012)}
DEC1 = {'Galileo I': (1, 2, 2, 1), 'NEAR': (1, 2, 2, 1), 'Cassini': (1, 2, 2, 1), 'Rosetta I': (2, 2, 2, 1),
        'Galileo II': (1, 2, 2, 1), 'MESSENGER': (2, 3, 3, 3)}
for r in table1(fl):
    k = r['spacecraft']
    for lab, key, i in (('dE', 'dE', 0), ('dLz', 'dLz_1e5', 1), ('sigma', 'sigma_1e5', 2), ('dLz/J', 'dLz_over_J_1e6', 3)):
        check(f"Table 1 {k} {lab}", r[key], T1[k][i], DEC1[k][i])

# Table 4: rp, Vp, eps, dcos
T4 = {'Galileo I': (7.331, 13.741, 2.473, 0.1486), 'NEAR': (6.910, 12.740, 1.814, 0.6254), 'Cassini': (7.546, 19.025, 5.852, -0.0215),
      'Rosetta I': (8.327, 10.519, 1.312, 0.1726), 'MESSENGER': (8.718, 10.387, 1.360, 0.0044), 'Rosetta II': (11.678, 12.487, 3.568, 0.0344),
      'EPOXI I': (21.945, 6.944, 1.655, 0.0391), 'EPOXI II': (49.828, 5.289, 2.497, 0.5197), 'Rosetta III': (8.859, 13.340, 2.955, 0.0375),
      'EPOXI III': (36.859, 5.727, 2.033, -0.5499), 'Juno': (6.933, 14.931, 2.877, 0.1977), 'Galileo II': (6.674, 14.080, 2.319, -0.1699)}
for r in table4(fl):
    k = r['spacecraft']
    for lab, key, i, dec in (('rp', 'rp_Mm', 0, 3), ('Vp', 'Vp_km_s', 1, 3), ('eps', 'eps', 2, 3), ('dcos', 'dcos', 3, 4)):
        check(f"Table 4 {k} {lab}", r[key], T4[k][i], dec)

# Table 5
T5 = {'Galileo I': (7.19, 0.55, 2), 'NEAR': (3.97, 0.003, 3), 'Cassini': (47.9, 24.0, 1), 'Rosetta I': (1.58, 0.03, 2),
      'Galileo II': (6.06, 1.46, 2), 'MESSENGER': (0.72, 0.41, 2)}
for r in single_encounter_couplings(fl):
    k = r['spacecraft']; lam, sig, dec = T5[k]
    check(f"Table 5 {k} lambda1 (1e3)", r['lam_1e3'], lam, 1 if k == 'Cassini' else 2)
    check(f"Table 5 {k} sigma (1e3)", r['sigma_1e3'], sig, dec)

# Eq. (16) and the quoted statistics
fit = main_fit(fl)
check('Eq. (16) lambda1 (1e3)', fit['lam'] / 1e3, 3.902, 3)
check('Eq. (16) standard error (1e3)', fit['se'] / 1e3, 0.576, 3)
check('Eq. (16) R^2 about the mean', fit['r2_mean'], 0.904, 3)
check('Eq. (16) R^2 about the origin', fit['r2_origin'], 0.939, 3)
check('Eq. (16) t', fit['t'], 6.77, 2)
check('Eq. (16) p', fit['p'], 0.0066, 4)
check('axial displacement d (m)', fit['d_m'], 17.3, 1)
check('sigma of d (m)', fit['sigma_d_m'], 2.6, 1)

# Table 6
T6 = {'Galileo I': (1.389, -2.45, 5.23), 'NEAR': (7.092, 2.25, 11.93), 'Cassini': (-0.137, -5.38, 5.11),
      'Rosetta I': (1.632, -1.37, 4.63), 'Galileo II': (-1.869, -5.85, 2.11)}
rows, t975 = predictions(fl, fit)
check('Table 6 t_0.025,3', t975, 3.18, 2)
for r in rows:
    k = r['spacecraft']
    check(f"Table 6 {k} ASTG", r['astg'], T6[k][0], 3)
    check(f"Table 6 {k} interval low", r['lo'], T6[k][1], 2)
    check(f"Table 6 {k} interval high", r['hi'], T6[k][2], 2)

# Table 7 and PRESS
T7 = {'Galileo I': (0.030, 1.00, 0.03, 3802), 'NEAR': (0.900, 0.34, 1.03, 3318), 'Cassini': (0.0002, -0.94, None, 3896),
      'Rosetta I': (0.070, -1.11, 0.09, 4077)}
loo, press, rss = leave_one_out(fl, fit)
for r in loo:
    k = r['encounter']; h, sr, cd, lo = T7[k]
    check(f"Table 7 {k} leverage", r['leverage'], h, 4 if k == 'Cassini' else 3)
    check(f"Table 7 {k} standardised residual", r['std_residual'], sr, 2)
    if cd is None:
        ok = r['cooks_d'] < 0.01; results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {'Table 7 Cassini Cook D below 0.01':55s} computed {r['cooks_d']:.5f}")
    else:
        check(f"Table 7 {k} Cook's D", r['cooks_d'], cd, 2)
    check(f"Table 7 {k} leave-one-out lambda1", r['lam_loo'], lo, 0)
check('PRESS (mm/s)^2', press, 6.3, 1)
check('in-sample RSS (mm/s)^2', rss, 4.7, 1)

# Eq. (17)
c = control_fit(fl)
check('Eq. (17) lambda1 (1e3)', c['lam'] / 1e3, 4.007, 3)
check('Eq. (17) standard error (1e3)', c['se'] / 1e3, 0.539, 3)
check('Eq. (17) R^2', c['r2_mean'], 0.922, 3)
check('Eq. (17) t', c['t'], 7.44, 2)
check('Eq. (17) p', c['p'], 0.0017, 4)

# Table 8 and chi^2
T8 = {'MESSENGER': (0.038, 0.022, 4.5, 3.6), 'Rosetta II': (0.139, 0.749, 1.4, 7.5), 'Rosetta III': (0.247, 0.766, 2.3, 7.7),
      'Juno': (1.901, 4.430, 5.5, 22.1), 'EPOXI I': (0.081, 0.207, 0.8, 2.1), 'EPOXI II': (0.273, 3.647, 2.5, 36.5),
      'EPOXI III': (-0.488, -3.326, 4.0, 33.3)}
rows8, lam, se = table8(fl)
for r in rows8:
    k = r['spacecraft']
    check(f"Table 8 {k} ASTG", r['astg'], T8[k][0], 3)
    check(f"Table 8 {k} Anderson", r['anderson'], T8[k][1], 3)
    check(f"Table 8 {k} ASTG deviation", r['z_astg'], T8[k][2], 1)
    check(f"Table 8 {k} Anderson deviation", r['z_anderson'], T8[k][3], 1)
check('chi^2 over the nulls, ASTG', sum(r['z_astg']**2 for r in rows8), 81, 0)
check('chi^2 over the nulls, Anderson', sum(r['z_anderson']**2 for r in rows8), 3059, 0)

n_fail = results.count(False)
print(f"\n{len(results) - n_fail} of {len(results)} checks pass")
sys.exit(1 if n_fail else 0)
