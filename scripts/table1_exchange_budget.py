"""Table 1 of the paper: energy and angular momentum exchanged with the Earth's rotation,
Delta E = Vinf * Delta Vinf and Delta L_z = Delta E / omega (Eq. 7), from the measured Delta Vinf and sigma
of Anderson et al. (2008).  J is the specific orbital angular momentum, J = sqrt(GM * rp (1 + eps))."""
import numpy as np
import astg

ORDER = ['Galileo I', 'NEAR', 'Cassini', 'Rosetta I', 'Galileo II', 'MESSENGER']


def table1(flybys=None):
    flybys = flybys or astg.load_flybys()
    rows = []
    for k in ORDER:
        f = flybys[k]; d = astg.derived(f)
        de = d['vinf'] * f['dVinf_mm_s'] * 1e-3                 # J/kg
        dl = de / astg.OMEGA                                     # m^2/s
        sig = d['vinf'] * f['sigma_dVinf_mm_s'] * 1e-3 / astg.OMEGA
        j = np.sqrt(astg.GM * d['semilatus'])
        rows.append(dict(spacecraft=k, dE=de, dLz_1e5=dl / 1e5, sigma_1e5=sig / 1e5, dLz_over_J_1e6=dl / j * 1e6))
    return rows


if __name__ == '__main__':
    print('Table 1. Exchange with the Earth\'s rotation (Anderson et al. 2008 Delta Vinf)')
    print(f"{'Spacecraft':12s} {'dE (J/kg)':>10s} {'dLz (1e5 m2/s)':>15s} {'sigma (1e5)':>12s} {'dLz/J (1e-6)':>13s}")
    for r in table1():
        print(f"{r['spacecraft']:12s} {r['dE']:+10.2f} {r['dLz_1e5']:+15.3f} {r['sigma_1e5']:12.3f} {r['dLz_over_J_1e6']:+13.3f}")
