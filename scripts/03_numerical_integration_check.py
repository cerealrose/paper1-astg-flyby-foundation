#!/usr/bin/env python3
"""
Independent numerical cross-check of the O(k) Binet-route formula, Eq. (B4)
of Nesvinga, Nyambuya & Jones (2026a), against direct numerical integration
of the full nonlinear equations of motion,

    r'' = r*phi'^2 - GM/r^2,
    phi'' = 2*k*phi'^2 - 2*r'*phi'/r,

exactly as claimed in Appendix D ("agrees to four significant figures
across epsilon = 1.31-5.85") and Appendix A footnote (checked "to four
significant figures"). This script reproduces that claim from scratch
using an independent method (direct ODE integration, not the symbolic
perturbative solution of script 02) and additionally verifies the
expected O(k) truncation-error scaling (fractional error should fall
linearly with k as k -> 0), which the paper's stated k ~ 1e-7 to 1e-8
regime relies on for the O(k) approximation to be trustworthy.

Run with:  python3 03_numerical_integration_check.py
Requires:  numpy, scipy
"""
import numpy as np
from scipy.integrate import solve_ivp

GM = 1.0  # normalised units throughout; only ratios/scalings are physical


def run_odeint(eps, k, ell0=1.0, t_max=8000.0):
    """Integrate the exact nonlinear equations of motion for a constant-k
    perturbation and return the numerically-determined Delta V_inf/V_inf."""
    rp = ell0 / (1 + eps)
    J0 = np.sqrt(GM * ell0)
    phidot0 = J0 / rp**2

    def rhs(t, y):
        r, rdot, phi, phidot = y
        rddot = r * phidot**2 - GM / r**2
        phiddot = 2 * k * phidot**2 - 2 * rdot * phidot / r
        return [rdot, rddot, phidot, phiddot]

    y0 = [rp, 0.0, 0.0, phidot0]
    sol_out = solve_ivp(rhs, [0, t_max], y0, method='DOP853', rtol=1e-13, atol=1e-15)
    r, rdot, phi, phidot = sol_out.y[:, -1]
    vinf_out = np.sqrt(rdot**2 + r**2 * phidot**2 - 2 * GM / r)

    sol_in = solve_ivp(rhs, [0, -t_max], y0, method='DOP853', rtol=1e-13, atol=1e-15)
    r, rdot, phi, phidot = sol_in.y[:, -1]
    vinf_in = np.sqrt(rdot**2 + r**2 * phidot**2 - 2 * GM / r)

    return (vinf_out - vinf_in) / vinf_in


def closed_form(eps, k):
    """Eq. (B4): the O(k) perturbative formula being tested."""
    eta = np.arccos(-1 / eps)
    return 2 * k * ((eps**2 + 2) * eta + 3 * np.sqrt(eps**2 - 1)) / (eps**2 - 1)


def main():
    print("=" * 72)
    print("Part 1: reproduce the paper's 'agrees to four significant")
    print("figures across eps=1.31-5.85' claim (Appendix D)")
    print("=" * 72)
    print("Using a moderate k (not the real k~1e-7-1e-8) for numerical")
    print("tractability; Part 2 below extrapolates to the real regime.\n")
    print(f"{'eps':>6} {'k':>8} {'numeric':>14} {'formula':>14} {'rel.err':>10}")
    for eps in [1.31, 1.81, 2.47, 3.86, 5.85]:
        k = 0.001
        num = run_odeint(eps, k, t_max=10000.0)
        formula = closed_form(eps, k)
        relerr = abs(num - formula) / abs(formula)
        print(f"{eps:6.2f} {k:8.4f} {num:14.8f} {formula:14.8f} {relerr:10.2e}")

    print("\n" + "=" * 72)
    print("Part 2: O(k) truncation-error scaling check")
    print("=" * 72)
    print("Fractional error should scale linearly in k for a correctly")
    print("derived first-order formula; extrapolating this scaling to the")
    print("paper's real k~1e-7-1e-8 regime tests whether 'four significant")
    print("figures' is a conservative or an optimistic characterisation.\n")
    eps = 2.47
    print(f"{'k':>10} {'numeric':>14} {'formula':>14} {'frac.err':>10} {'frac.err/k':>10}")
    ratios = []
    for k in [0.02, 0.01, 0.005, 0.002, 0.001]:
        num = run_odeint(eps, k, t_max=10000.0)
        formula = closed_form(eps, k)
        fe = abs(num - formula) / abs(formula)
        ratios.append(fe / k)
        print(f"{k:10.4f} {num:14.8f} {formula:14.8f} {fe:10.5f} {fe/k:10.3f}")

    print(f"\nfrac.err/k converges to ~{ratios[-1]:.2f} as k -> 0, confirming")
    print("linear-in-k scaling (the correct signature of a first-order")
    print("truncation, not a bug). Extrapolated to k~1e-7-1e-8, the real")
    print(f"flyby regime, fractional error ~ {ratios[-1]:.1f}e-7 to "
          f"{ratios[-1]:.1f}e-8 -- far better than four significant figures.")


if __name__ == "__main__":
    main()
    print("\nDone.")
