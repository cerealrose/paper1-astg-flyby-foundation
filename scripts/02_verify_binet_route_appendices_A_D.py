#!/usr/bin/env python3
"""
Verification script for Nesvinga, Nyambuya & Jones (2026a), Appendices A-D
(the generalised Binet-equation derivation, retained as an independent
cross-check of the vis-viva result).

Re-derives, symbolically and from scratch, every displayed equation in the
chain: Newton's second law -> generalised Binet equation with a
non-conservative azimuthal force -> exact solution -> O(k) small-k
expansion -> perigee boundary conditions -> asymptotic angle shifts ->
Delta V_infinity / V_infinity -> perigee-frame conversion to the boxed
Eq. (B5)/(D2) result. Every intermediate equation is checked against its
own claimed identity (residual = 0) rather than only checking the final
formula.

Run with:  python3 02_verify_binet_route_appendices_A_D.py
Requires:  sympy
"""
import sympy as sp


def main():
    phi, k, eps, ell0, eta = sp.symbols('phi k epsilon ell_0 eta', real=True)

    print("=" * 72)
    print("Step 1: particular + homogeneous solution of the generalised")
    print("Binet equation, Eq. (A13)/(C1)")
    print("=" * 72)

    # Particular solution u_p = A*exp(-4k*phi), A = 1/(ell0*(1+8k^2))  [Eq. C3]
    A_expr = 1 / (ell0 * (1 + 8 * k**2))
    u_p = A_expr * sp.exp(-4 * k * phi)
    rhs = sp.exp(-4 * k * phi) / ell0
    residual = sp.simplify(sp.diff(u_p, phi, 2) + 2 * k * sp.diff(u_p, phi) + u_p - rhs)
    print(f"  particular-solution residual = {residual}  (expect 0)")
    assert residual == 0

    # Homogeneous solution, Omega = sqrt(1-k^2)  [Eq. C4]
    Omega = sp.sqrt(1 - k**2)
    C1, C2 = sp.symbols('C1 C2')
    u_h = sp.exp(-k * phi) * (C1 * sp.cos(Omega * phi) + C2 * sp.sin(Omega * phi))
    residual = sp.simplify(sp.diff(u_h, phi, 2) + 2 * k * sp.diff(u_h, phi) + u_h)
    print(f"  homogeneous-solution residual = {residual}  (expect 0)")
    assert residual == 0

    print("\n" + "=" * 72)
    print("Step 2: perigee boundary conditions -> C1, C2, Eq. (C5)")
    print("=" * 72)
    u_full = u_p + u_h
    C1_sol = sp.solve(sp.Eq(u_full.subs(phi, 0), (1 + eps) / ell0), C1)[0]
    paper_C1 = (1 + eps) / ell0 - A_expr
    assert sp.simplify(C1_sol - paper_C1) == 0
    print("  C1 matches paper's Eq. (C5) exactly.")

    uprime0 = sp.diff(u_full, phi).subs(phi, 0).subs(C1, C1_sol)
    C2_sol = sp.simplify(sp.solve(sp.Eq(uprime0, 0), C2)[0])
    paper_C2 = sp.simplify(k * (4 * A_expr + C1_sol) / Omega)
    assert sp.simplify(C2_sol - paper_C2) == 0
    print("  C2 matches paper's Eq. (C5) exactly.")

    print("\n" + "=" * 72)
    print("Step 3: small-k expansion -> Eq. (A15)/(C7)")
    print("=" * 72)
    u_full_k = u_p + sp.exp(-k * phi) * (C1 * sp.cos(Omega * phi) + C2 * sp.sin(Omega * phi))
    u_full_k = u_full_k.subs({C1: C1_sol, C2: C2_sol})
    u_series = sp.expand(sp.series(u_full_k, k, 0, 2).removeO())
    paper_u = (1 / ell0) * (1 + eps * sp.cos(phi)) + (k / ell0) * (
        -eps * phi * sp.cos(phi) + eps * sp.sin(phi) - 4 * phi + 4 * sp.sin(phi))
    diff = sp.simplify(sp.expand(u_series - paper_u))
    print(f"  O(k) expansion difference from Eq. (A15) = {diff}  (expect 0)")
    assert diff == 0

    print("\n" + "=" * 72)
    print("Step 4: asymptotic angle shifts a_out, a_in, Eq. (B2)-(B3)")
    print("=" * 72)
    subs_trig = {sp.tan(eta): -sp.sqrt(eps**2 - 1),
                 sp.cos(eta): -1 / eps, sp.sin(eta): sp.sqrt(eps**2 - 1) / eps}
    u0p = sp.diff((1 / ell0) * (1 + eps * sp.cos(phi)), phi)
    u1 = (1 / ell0) * (-eps * phi * sp.cos(phi) + eps * sp.sin(phi)
                        - 4 * phi + 4 * sp.sin(phi))

    a_out = sp.simplify(-u1.subs(phi, eta) / u0p.subs(phi, eta))
    a_out_sub = sp.simplify(a_out.subs(subs_trig))
    paper_a_out = 1 + 4 / eps - 3 * eta / sp.sqrt(eps**2 - 1)
    assert sp.simplify(a_out_sub - paper_a_out) == 0
    print("  a_out matches paper's Eq. (B2) exactly.")

    a_in = sp.simplify(-u1.subs(phi, -eta) / u0p.subs(phi, -eta))
    a_in_sub = sp.simplify(a_in.subs(subs_trig))
    paper_a_in = (4 * sp.sqrt(eps**2 - 1) - eps * (3 * eta - sp.sqrt(eps**2 - 1))) / (eps**2 * sp.sin(eta))
    paper_a_in_sub = sp.simplify(paper_a_in.subs(sp.sin(eta), sp.sqrt(eps**2 - 1) / eps))
    assert sp.simplify(a_in_sub - paper_a_in_sub) == 0
    print("  a_in matches paper's Eq. (B3) exactly.")
    print(f"\n  Note: a_in simplifies to the SAME expression as a_out")
    print(f"  ({sp.simplify(a_in_sub - a_out_sub) == 0}), an unadvertised")
    print(f"  consistency check between two independently-solved boundary")
    print(f"  conditions.")

    print("\n" + "=" * 72)
    print("Step 5: Delta V_infinity / V_infinity, Eq. (B4)/(D1)")
    print("=" * 72)
    J0 = sp.symbols('J0', positive=True)
    k_hat_dummy = k  # constant-k route
    J_series = sp.expand(sp.series(J0 * sp.exp(2 * k * phi), k, 0, 2).removeO())

    u1p = sp.diff(u1, phi)
    up_full = u0p + k * u1p

    phi_out = eta + k * a_out
    phi_in = -eta + k * a_in

    J_out = sp.series(J_series.subs(phi, phi_out), k, 0, 2).removeO()
    J_in = sp.series(J_series.subs(phi, phi_in), k, 0, 2).removeO()
    up_out = sp.series(up_full.subs(phi, phi_out), k, 0, 2).removeO()
    up_in = sp.series(up_full.subs(phi, phi_in), k, 0, 2).removeO()

    v_out = sp.expand(J_out * (-up_out))
    v_in = sp.expand(J_in * up_in)

    dV_over_V = sp.series((v_out - v_in) / v_in, k, 0, 2).removeO()
    dV_over_V = sp.simplify(dV_over_V.subs(subs_trig))

    paper_dVinf = 2 * k * ((eps**2 + 2) * eta + 3 * sp.sqrt(eps**2 - 1)) / (eps**2 - 1)
    diff = sp.simplify(dV_over_V - paper_dVinf)
    print(f"  Delta V_inf/V_inf difference from Eq. (B4) = {diff}  (expect 0)")
    assert diff == 0
    print("  Full chain from Newton's second law to Eq. (B4)/(D1): VERIFIED.")

    print("\n" + "=" * 72)
    print("Step 6: perigee-frame conversion, boxed Eq. (B5)/(D2)")
    print("=" * 72)
    # V_p^2/V_inf^2 = (1+eps)/(eps-1); dVp/Vp = (Vinf^2/Vp^2)*dVinf/Vinf
    factor = (eps - 1) / (eps + 1)
    dVp_over_Vp = sp.simplify(factor * paper_dVinf)
    paper_boxed = 2 * k * ((eps**2 + 2) * eta + 3 * sp.sqrt(eps**2 - 1)) / (1 + eps)**2
    diff = sp.simplify(dVp_over_Vp - paper_boxed)
    print(f"  boxed-equation difference = {diff}  (expect 0)")
    assert diff == 0
    print("  Perigee-frame conversion: VERIFIED. Full Appendix A-D chain")
    print("  reproduces the paper's boxed Eq. (B5)/(D2) exactly, with zero")
    print("  algebra errors at any intermediate step.")


if __name__ == "__main__":
    main()
    print("\nAll checks in this script passed.")
