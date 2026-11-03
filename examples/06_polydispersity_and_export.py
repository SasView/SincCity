"""Size polydispersity, and the ways a curve leaves the calculator.

The core-shell particle of core_shell_sphere.sld in D2O, with a Gaussian
spread of sizes (sasmodels' radius_pd for the whole structure). A spread
fills in the minima of the form factor. Each curve is then written as CSV,
converted to a sasdata Data1D for a fit (when sasdata is installed), and
plotted together (with --plot, when plotly is installed).

    python examples/06_polydispersity_and_export.py [--plot]
"""

import sys
from pathlib import Path

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'
OUTPUT = Path(__file__).resolve().parent / 'output'

sinccity.set_verbosity('quiet')
gsc = GenericScatteringCalculator(DATA / 'core_shell_sphere.sld')
# The grid stops at 0.25 1/Ang: the largest size factor (1 + 3 size_pd) times
# q_max times the 3 Ang pixels stays below 1, so treating pixels as points holds.
gsc.set_q_grid(0.01, 0.25, 120, spacing='linear')
gsc.set_params(sld_solvent=6.36, scale=0.01)  # a 1 % volume fraction, no background

# size_pd is the relative standard deviation of the size; the one spectrum is
# evaluated at f q for every size factor f, so a spread costs no pair pass.
spreads = (0.0, 0.03, 0.06, 0.1)
profiles = gsc.calculate_series(
    [{'label': f'size_pd {pd:g}', 'params': {'size_pd': pd}} for pd in spreads]
)

# The first minimum of the monodisperse curve, and how deep each curve is there.
ideal = profiles[0].intensity
first = 1 + int(np.argmax((ideal[1:-1] < ideal[:-2]) & (ideal[1:-1] < ideal[2:])))
q_min = gsc.q[first]
print(f'first minimum of the monodisperse curve at q = {q_min:.4f} 1/Ang')
print(f'{"curve":<14} {"sizes":>6} {"I(q_min) / I(q_0)":>18} {"time (s)":>9}')
for profile in profiles:
    depth = profile.intensity[first] / profile.intensity[0]
    print(f'{profile.label:<14} {profile.n_sizes:>6} {depth:18.3e} {profile.elapsed:9.3f}')
depths = [profile.intensity[first] / profile.intensity[0] for profile in profiles]
assert all(a < b for a, b in zip(depths, depths[1:])), 'a wider spread should fill the minimum more'
assert profiles[0].n_sizes == 1 and all(profile.n_sizes > 1 for profile in profiles[1:])

# CSV: q and I(q) with the provenance in an ASCII header.
OUTPUT.mkdir(exist_ok=True)
for profile in profiles:
    path = OUTPUT / f'core_shell_{profile.label.replace(" ", "_")}.csv'
    profile.save_csv(path)
print(f'\nsaved {len(profiles)} CSV files in {OUTPUT}; header of the last one:')
print(path.read_text().splitlines()[0][:160] + ' ...')
saved = np.loadtxt(path, delimiter=',')  # the header lines are comments
assert np.allclose(saved, np.column_stack([profiles[-1].q, profiles[-1].intensity])), 'CSV'

# A sasdata Data1D, for IQuana.set_data or any sasmodels fit. The curve is
# theory, so uncertainties are given explicitly, here 2 %.
try:
    data = profiles[2].to_data1d(dI=0.02 * profiles[2].intensity)
except ImportError as exc:
    print(f'\nno Data1D: {exc}')
else:
    print(f'\nData1D "{data.title}": {data.x.size} points, q {data.qmin:g} to {data.qmax:g}')
    print(f'process: {data.process[-1].description[:100]} ...')

if '--plot' in sys.argv:
    fig = profiles[0].plot(show=False)
    for profile in profiles[1:]:
        fig.add_scatter(x=profile.q, y=profile.intensity, mode='lines', name=profile.label)
    fig.data[0].name = profiles[0].label
    fig.show()
