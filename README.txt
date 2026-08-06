================================================================================
 Plausible Origins of the Earth Flyby Anomalies at the Perigee Point
 P. A. Nesvinga, G. G. Nyambuya, B. C. Jones
 NUST Fundamental Theoretical Physics Group
================================================================================

This repository accompanies the manuscript "Plausible Origins of the Earth
Flyby Anomalies at the Perigee Point" (submitted to Physical Review D). It
contains the manuscript source and a complete, independently re-run set of
verification scripts covering every derivation this paper claims to be
"verified" or "checked" in its text -- so a referee can reproduce every
such claim directly, rather than take it on faith.

This is the theory-foundation paper of a companion series. The full
statistical confrontation of the formula derived here with the complete
Earth-flyby dataset is carried out in the companion paper "Earth-Flyby
Velocity Anomalies within the Azimuthally Symmetric Theory of Gravitation"
(separate repository). Dedicated treatments of the 2013 Juno flyby and the
December 2026 Europa Clipper prediction are likewise separate, forthcoming
companion papers, not included here.

--------------------------------------------------------------------------------
REPOSITORY CONTENTS
--------------------------------------------------------------------------------

paper/
    paper_1_corrected.tex   Full manuscript source (REVTeX 4-2, PRD format)
    thebibliography.bib     Bibliography (BibTeX)

scripts/
    01_verify_visviva_derivation.py
        Symbolically re-derives Eqs. (10)-(13), the paper's central result,
        from the stated ell=1 ASTG potential via the exact vis-viva relation.
        Also verifies, via direct substitution into the angular Laplace
        equation, the claim underlying Sec. IX.A: that cos(theta) is the
        unique ell=1 solution under axisymmetry, and that cos(delta) (the
        empirically-preferred variable) solves none of ell=0..3 -- the
        basis for the paper stating that the Delta cos(delta) angular form
        is adopted on empirical/structural grounds, not derived from the
        stated potential. Finally reproduces the R^2=0.90 vs 0.53 fit
        comparison quoted in the same section.

    02_verify_binet_route_appendices_A_D.py
        Full symbolic re-derivation, from Newton's second law, of the
        entire generalised-Binet-equation cross-check (Appendices A-D):
        the non-conservative azimuthal-force equation of motion, its exact
        solution, the O(k) small-k expansion, the perigee boundary
        conditions, the asymptotic angle shifts, and the final boxed
        result Eq. (B5)/(D2). Every intermediate equation is checked
        against its own claimed identity, not only the final formula.

    03_numerical_integration_check.py
        Independent numerical cross-check of the O(k) formula via direct
        integration (SciPy, DOP853) of the full nonlinear equations of
        motion -- reproducing the paper's claim (Appendix D) that the
        perturbative formula "agrees to four significant figures across
        epsilon=1.31-5.85", and additionally verifying the expected
        linear-in-k truncation-error scaling, which is what justifies
        trusting the O(k) approximation at the real flyby k~1e-7-1e-8.

    04_verify_appendix_E_phidep.py
        Symbolic derivation and numerical cross-check of Appendix E's
        exact O(k_hat) integration of the non-autonomous (phi-dependent
        coupling) Binet equation, Eq. (E5) -- the most delicate result in
        the paper. Includes both the symbolic derivation of u_1(phi),
        Eq. (E4), with boundary conditions verified exactly, and a direct
        numerical ODE integration cross-check with a truncation-error
        scaling test, exactly as used to validate Eq. (E5) before it was
        included in the manuscript.

    05_lambda2_suppression_estimate.py
        Reproduces the order-of-magnitude lambda_2/lambda_1 estimate used
        in Sec. XI's discussion of the Juno tension, using Nyambuya
        (2015)'s phenomenological scaling law applied to Earth's own
        (correctly-defined) specific spin angular momentum. Includes an
        explicit robustness check showing the paper's conclusion (the l=2
        term is ~10 orders of magnitude too small to resolve Juno) holds
        under either candidate definition of Earth's specific angular
        momentum.

requirements.txt   Python dependencies (numpy, scipy, sympy)
LICENSE            MIT License (applies to all code in scripts/)
CITATION.cff       Machine-readable citation metadata

--------------------------------------------------------------------------------
REPRODUCING THE RESULTS
--------------------------------------------------------------------------------

1. Install dependencies:
       pip install -r requirements.txt

2. Run any script directly, e.g.:
       python3 scripts/01_verify_visviva_derivation.py

   Each script is self-contained, prints its own step-by-step derivation
   and pass/fail status for every assertion, and exits with a non-zero
   status code if any check fails. Scripts 01-02 run in seconds; 03-04
   involve numerical ODE integration and run in under two minutes each on
   a standard laptop.

3. To compile the manuscript itself (requires a LaTeX distribution with
   REVTeX 4-2, e.g. TeX Live's texlive-publishers package):
       cd paper
       pdflatex paper_1_corrected.tex
       bibtex paper_1_corrected
       pdflatex paper_1_corrected.tex
       pdflatex paper_1_corrected.tex

--------------------------------------------------------------------------------
WHAT THIS REPOSITORY DOES NOT CLAIM
--------------------------------------------------------------------------------

In the interest of transparency to referees: this repository verifies that
the manuscript's derivations are internally correct and that its numerical
claims are reproducible. It does not resolve two open questions the
manuscript itself states plainly (Sec. IX.A and Sec. XI):

  (i)  The empirically-adopted Delta cos(delta) angular dependence is not
       derived from the stated ell=1 potential -- script 01 verifies this
       is a structural fact (via the Legendre-equation check), not an
       oversight to be fixed by more careful algebra.

  (ii) The l=2 suppression estimate in script 05 is an order-of-magnitude
       plausibility argument, not a rigorous bound -- it extrapolates a
       scaling law fit to Solar perihelion precession to Earth's flyby
       regime, an extrapolation the manuscript explicitly does not claim
       is justified.

--------------------------------------------------------------------------------
CONTACT
--------------------------------------------------------------------------------
P. A. Nesvinga <praise.nesvinga@nust.ac.zw>
Department of Applied Physics, National University of Science & Technology,
Bulawayo, Zimbabwe
