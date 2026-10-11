"""Table 4 of the paper: derived orbital quantities from Eq. (9),
rp = R_earth + hp, Vp = sqrt(Vinf^2 + 2 GM / rp), eps = 1 + rp Vinf^2 / GM, and dcos."""
import astg

ORDER = ['Galileo I', 'NEAR', 'Cassini', 'Rosetta I', 'MESSENGER', 'Rosetta II', 'EPOXI I', 'EPOXI II',
         'Rosetta III', 'EPOXI III', 'Juno', 'Galileo II']


def table4(flybys=None):
    flybys = flybys or astg.load_flybys()
    rows = []
    for k in ORDER:
        d = astg.derived(flybys[k])
        rows.append(dict(spacecraft=k, rp_Mm=d['rp'] / 1e6, Vp_km_s=d['vp'] / 1e3, eps=d['eps'], dcos=d['dcos']))
    return rows


if __name__ == '__main__':
    print('Table 4. Derived orbital quantities')
    print(f"{'Spacecraft':12s} {'rp (Mm)':>8s} {'Vp (km/s)':>10s} {'eps':>7s} {'dcos':>8s}")
    for r in table4():
        print(f"{r['spacecraft']:12s} {r['rp_Mm']:8.3f} {r['Vp_km_s']:10.3f} {r['eps']:7.3f} {r['dcos']:+8.4f}")
