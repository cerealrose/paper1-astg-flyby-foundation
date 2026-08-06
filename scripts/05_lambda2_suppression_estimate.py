#!/usr/bin/env python3
"""
Reproduces the order-of-magnitude lambda_2/lambda_1 estimate quoted in
Nesvinga, Nyambuya & Jones (2026a), Sec. XI ("Future Direction"), item 3
("The Juno tension").

Uses Nyambuya (2015)'s phenomenological scaling law,
    lambda_{1,l} = (-1)^(l+1) * (S_star / S*_l)^l,
applied to Earth's own specific spin angular momentum, to estimate the
relative size of the l=2 term compared to l=1, and combines this with the
intrinsic (GM/(r_p c0^2)) suppression of the l=2 multipole term to bound
the l=2 contribution to the Juno discrepancy.

An earlier, less careful pass at this calculation (not by the authors)
silently used two different, inconsistent definitions of Earth's specific
angular momentum -- one dimensionally inconsistent (giving a nonsensical
result), one an unjustified idealised-shell approximation. This script
uses the single, physically correct definition (Earth's actual measured
moment-of-inertia factor) and shows the paper's stated conclusion is
insensitive to the choice in any case.

Run with:  python3 05_lambda2_suppression_estimate.py
Requires:  numpy
"""
import numpy as np


def main():
    # Physical constants and Earth parameters (CODATA / standard values)
    M_E = 5.972e24       # kg
    R_E = 6.371e6         # m
    omega_E = 7.292115e-5  # rad/s (sidereal), same value used throughout
                            # this paper series, sourced from Anderson (2008)
    G = 6.67430e-11       # m^3 kg^-1 s^-2
    c0 = 2.99792458e8     # m/s

    # Nyambuya (2015) fitted values, Eq. (7) of Paper I:
    # G11 = (678.00+/-4.00)*G, G12 = (6.06+/-0.05)*G
    G11 = 678.0 * G
    G12 = 6.06 * G

    print("=" * 72)
    print("Step 1: dimensional check of S*_l = sqrt(G_l * M * R)")
    print("=" * 72)
    S_star_1 = np.sqrt(G11 * M_E * R_E)
    S_star_2 = np.sqrt(G12 * M_E * R_E)
    print(f"  S*_1 = {S_star_1:.4e}  [units: m^2/s, i.e. SPECIFIC ang. mom.]")
    print(f"  S*_2 = {S_star_2:.4e}")
    print("  For (S_star/S*_l) to be dimensionless, S_star must ALSO be")
    print("  specific angular momentum (m^2/s), not total (kg m^2/s).")

    print("\n" + "=" * 72)
    print("Step 2: two candidate definitions of Earth's S_star")
    print("=" * 72)
    I_E = 0.3307 * M_E * R_E**2  # Earth's actual, measured moment of inertia
    S_total = I_E * omega_E
    print(f"  Total angular momentum S=I*omega = {S_total:.4e} kg m^2/s")
    print("  (dimensionally WRONG for this formula -- would give a")
    print(f"   nonsensical lambda_11 ~ {S_total/S_star_1:.2e}; discarded)")

    S_idealised = R_E**2 * omega_E  # idealised shell, k=1 (NOT Earth's real I/M)
    S_correct = 0.3307 * R_E**2 * omega_E  # physically correct: (I/M)*omega

    print(f"\n  Idealised-shell specific S = R_E^2*omega_E = {S_idealised:.4e} m^2/s")
    print("  (uses k=1, NOT Earth's actual moment-of-inertia factor --")
    print("   not used in the paper's final estimate)")
    print(f"\n  Physically correct specific S = (I_E/M_E)*omega_E")
    print(f"                                = 0.3307*R_E^2*omega_E")
    print(f"                                = {S_correct:.4e} m^2/s")
    print("  (this is the value used in the paper: Earth's REAL, measured")
    print("   moment-of-inertia factor, not an idealised approximation)")

    print("\n" + "=" * 72)
    print("Step 3: lambda_12/lambda_11 ratio for both candidate definitions")
    print("=" * 72)
    for label, S in [("idealised shell (k=1)", S_idealised),
                      ("physically correct (k=0.3307)", S_correct)]:
        lam11 = S / S_star_1
        lam12 = -(S / S_star_2)**2
        ratio = lam12 / lam11
        print(f"  [{label}]")
        print(f"    lambda_11 = {lam11:.4e}, lambda_12 = {lam12:.4e}, "
              f"ratio = {ratio:.4f}")

    ratio_used = -(S_correct / S_star_2)**2 / (S_correct / S_star_1)
    print(f"\n  Paper uses the physically correct value: "
          f"lambda_12/lambda_11 = {ratio_used:.3f}")
    print(f"  (rounds to the paper's stated '~-0.08')")

    print("\n" + "=" * 72)
    print("Step 4: combine with the intrinsic l=2 suppression factor")
    print("=" * 72)
    GM = G * M_E
    rp_typical = 7.0e6  # m, typical Earth-flyby perigee radius
    suppression = GM / (rp_typical * c0**2)
    print(f"  GM_Earth/(r_p*c0^2) at r_p~7000 km = {suppression:.3e}")

    combined = abs(ratio_used) * suppression
    print(f"\n  Combined l=2/l=1 term suppression = "
          f"|lambda_12/lambda_11| * GM/(r_p c0^2)")
    print(f"                                    = {combined:.3e}")

    dV1 = 2.0  # mm/s, order of the l=1 (Juno) prediction
    dV2 = dV1 * combined
    print(f"\n  If the l=1 term gives ~{dV1} mm/s (Juno's predicted anomaly),")
    print(f"  the l=2 term contributes ~{dV2:.2e} mm/s -- "
          f"{np.log10(dV1/dV2):.0f} orders of")
    print(f"  magnitude below the observed discrepancy, confirming the")
    print(f"  paper's statement that a straightforward l=2 term is very")
    print(f"  unlikely to resolve Juno on its own.")

    print("\n" + "=" * 72)
    print("Robustness check: does the choice of S_star definition matter?")
    print("=" * 72)
    for label, S in [("idealised shell", S_idealised), ("physically correct", S_correct)]:
        ratio = -(S / S_star_2)**2 / (S / S_star_1)
        combined = abs(ratio) * suppression
        print(f"  [{label}]: combined suppression = {combined:.2e} "
              f"({np.log10(dV1/(dV1*combined)):.0f} orders of magnitude)")
    print("\n  Conclusion is robust to which definition is used: both give")
    print("  a suppression far too small for the l=2 term to resolve Juno.")


if __name__ == "__main__":
    main()
    print("\nDone.")
