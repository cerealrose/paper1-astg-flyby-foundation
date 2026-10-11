"""Figures 1 and 2 of the paper, written to ../figures/fig1.png and ../figures/fig2_horizons.png.
Figure 1: diagnostics of the global fit (regression through the origin, predicted against observed,
standardised residuals, leave-one-out refits).  Figure 2: deviation of the seven null encounters from the
ASTG law and from the Anderson relation."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import astg
from null_encounters import table8

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'figures')
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 11,
                     'axes.spines.top': False, 'axes.spines.right': False, 'mathtext.fontset': 'dejavusans'})
BLUE, RED = '#1f4e79', '#b22222'
LAB = {'Galileo I': 'Gal I', 'NEAR': 'NEAR', 'Cassini': 'Cas', 'Rosetta I': 'Ros I'}


def figure1(fl):
    an = astg.ANOMALOUS
    d = [astg.derived(fl[k]) for k in an]
    X = np.array([q['X'] for q in d]); Vp = np.array([q['vp'] for q in d])
    Y = np.array([fl[k]['dVp_mm_s'] * 1e-3 for k in an]) / Vp
    fit = astg.fit_through_origin(X, Y); lam, se, s2 = fit['lam'], fit['se'], fit['s2']
    h = fit['leverage']; sres = fit['res'] / np.sqrt(s2 * (1 - h))
    fig, ax = plt.subplots(2, 2, figsize=(3679 / 400, 2960 / 400), dpi=400)
    a = ax[0, 0]; a.axhline(0, color='0.8', lw=0.8); a.axvline(0, color='0.8', lw=0.8)
    xx = np.linspace(0, 1.6e-10, 50); a.plot(xx * 1e10, lam * xx * 1e7, 'k-', lw=1.5, label=r'$\hat\lambda_1=%.0f$' % lam)
    a.scatter(X * 1e10, Y * 1e7, color=BLUE, zorder=3)
    for k, x, y in zip(an, X, Y):
        a.annotate(LAB[k], (x * 1e10, y * 1e7), xytext=(6, -10), textcoords='offset points', fontsize=9)
    a.set_xlabel(r'$X_i=R_g\,\Delta\!\cos\delta_i/[2r_{p,i}(1+\varepsilon_i)]$  ($\times10^{-10}$)')
    a.set_ylabel(r'$Y_i=\Delta V_p/V_p$  ($\times10^{-7}$)'); a.set_title('(a) Regression through the origin', fontsize=11)
    a.legend(frameon=False, loc='lower right', fontsize=9)
    a = ax[0, 1]; pr = lam * X * Vp * 1e3; ob = Y * Vp * 1e3
    a.plot([-3, 9], [-3, 9], 'k--', lw=0.9, label='1:1'); a.scatter(ob, pr, color=BLUE, zorder=3)
    for k, x, y in zip(an, ob, pr):
        a.annotate(LAB[k], (x, y), xytext=(6, -10), textcoords='offset points', fontsize=9)
    a.set_xlim(-3, 9); a.set_ylim(-3, 9); a.set_xlabel(r'observed $\Delta V_p$ (mm$\,$s$^{-1}$)'); a.set_ylabel(r'predicted $\Delta V_p$ (mm$\,$s$^{-1}$)')
    a.set_title('(b) Predicted against observed', fontsize=11); a.legend(frameon=False, loc='upper left', fontsize=9)
    a = ax[1, 0]; fitted = lam * X
    a.axhline(0, color='k', lw=0.9)
    for yv in (2, -2):
        a.axhline(yv, color='0.6', ls=':', lw=0.9)
    a.scatter(fitted * 1e7, sres, color=BLUE, zorder=3)
    for k, x, y in zip(an, fitted, sres):
        a.annotate(LAB[k], (x * 1e7, y), xytext=(6, -10), textcoords='offset points', fontsize=9)
    a.set_ylim(-2.6, 2.6); a.set_xlabel(r'fitted value ($\times10^{-7}$)'); a.set_ylabel('standardised residual'); a.set_title('(c) Residuals against fitted value', fontsize=11)
    a = ax[1, 1]
    loo = [X[np.arange(4) != i] @ Y[np.arange(4) != i] / (X[np.arange(4) != i] @ X[np.arange(4) != i]) for i in range(4)]
    a.axhspan(lam - se, lam + se, color=RED, alpha=0.12, lw=0); a.axhline(lam, color=RED, lw=1.2, label=r'full sample $\hat\lambda_1=%.0f$' % lam)
    a.bar([LAB[k] for k in an], loo, color='0.7', edgecolor='0.3', width=0.6)
    a.set_ylabel(r'$\hat\lambda_1^{(-i)}$'); a.set_ylim(0, 4600); a.set_title('(d) Leave-one-out sensitivity', fontsize=11); a.legend(frameon=False, loc='upper left', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig1.png'), dpi=400); plt.close(fig)


def figure2(fl):
    rows, lam, se = table8(fl)
    names = [r['spacecraft'] for r in rows]
    za = [r['z_astg'] for r in rows]; zA = [r['z_anderson'] for r in rows]
    fig, a = plt.subplots(figsize=(3679 / 400, 2000 / 400), dpi=400)
    yv = np.arange(len(names))[::-1]
    a.plot(zA, yv + 0.15, 'o', mfc='white', mec='0.3', label='Anderson relation, Eq. (1)')
    a.plot(za, yv - 0.15, 'o', color=BLUE, label='ASTG law, Eq. (12)')
    a.axvline(3, color='0.5', ls='--', lw=0.9)
    a.set_xscale('log'); a.set_xlim(0.4, 100); a.set_yticks(yv); a.set_yticklabels(names)
    a.set_xlabel('deviation from observation (combined standard deviations)')
    a.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3, 5, 10, 20, 50, 100])); a.xaxis.set_major_formatter(FixedFormatter(['0.5', '1', '2', '3', '5', '10', '20', '50', '100']))
    a.legend(frameon=False, loc='upper right', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig2_horizons.png'), dpi=400); plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    fl = astg.load_flybys()
    figure1(fl); figure2(fl)
    print('wrote', os.path.normpath(os.path.join(OUT, 'fig1.png')), 'and fig2_horizons.png')
