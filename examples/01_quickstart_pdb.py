"""Quick start: I(q) of buckminsterfullerene (C60) from a PDB file, in vacuum.

Reads the PDB file, prints what was read, computes I(q) with the default
kernel, checks it against the exact pair sum and against the closed form of
a thin spherical shell, and saves the curve as CSV.

    python examples/01_quickstart_pdb.py           # tables only
    python examples/01_quickstart_pdb.py --plot    # also the 3D points and the curve (plotly)
"""

import sys
from pathlib import Path

import numpy as np

from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'
OUTPUT = Path(__file__).resolve().parent / 'output'

# c60.pdb has no element column and the residue name '???'. Without help the
# elements would be guessed from the atom names, with a warning; elements=
# states them per residue name, in atom order.
gsc = GenericScatteringCalculator()
structure = gsc.load_structure(DATA / 'c60.pdb', elements={'???': ['C'] * 60})
structure.describe()

# In vacuum (sld_solvent = 0, the default) every atom scatters with its own
# scattering length. The grid runs past the first minimum of the shell.
gsc.set_q_grid(0.01, 2.0, 60, spacing='linear')
profile = gsc.calculate()  # kernel 'auto': the distance histogram for atoms
profile.describe()

# The histogram kernel rounds distances to a small bin; the pair sum is exact.
# Compare on the scale of the curve, I(0): near a zero of I(q) any relative
# difference is large.
gsc.set_params(kernel='pairs')
exact = gsc.calculate()
deviation = np.max(np.abs(profile.intensity - exact.intensity)) / exact.intensity[0]
print(f'\nhistogram vs exact pair sum: largest difference {deviation:.1e} of I(0)')
assert deviation < 1e-4, 'the histogram kernel is documented as exact to about 1e-4'

# All 60 atoms sit on one sphere of radius R. A continuous shell of the same
# total scattering length gives I(q) = I(0) sinc^2(q R); the discrete atoms
# follow it until q is comparable with the inverse of their spacing.
radius = np.linalg.norm(structure.xyz - structure.xyz.mean(axis=0), axis=1).mean()
shell = exact.intensity[0] * np.sinc(exact.q * radius / np.pi) ** 2
shell /= np.sinc(exact.q[0] * radius / np.pi) ** 2
print(f'\nshell radius R = {radius:.3f} Ang; first zero of sinc(qR) at q = {np.pi / radius:.3f}')
print(f'{"q (1/Ang)":>10} {"I(q) (1/cm)":>13} {"shell":>13} {"ratio":>7}')
for i in range(0, exact.q.size, 8):
    q, value, model = exact.q[i], exact.intensity[i], shell[i]
    print(f'{q:10.3f} {value:13.5g} {model:13.5g} {value / model:7.3f}')
# Up to the first zero the atoms are indistinguishable from the shell.
inner = exact.q < np.pi / radius
shell_error = np.max(np.abs(exact.intensity[inner] - shell[inner])) / exact.intensity[0]
assert shell_error < 1e-3, f'C60 departs from a spherical shell by {shell_error:.1e} of I(0)'

OUTPUT.mkdir(exist_ok=True)
profile.save_csv(OUTPUT / 'c60.csv')  # two columns, the provenance in the header
print(f'\nsaved {OUTPUT / "c60.csv"}')

if '--plot' in sys.argv:
    gsc.plot_points()  # SasView's Draw Points, with the CONECT bonds
    profile.plot()
