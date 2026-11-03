"""A surfactant micelle from a molecular-dynamics frame, in four contrasts.

micelle.gro is a frame as GROMACS writes it: 24 cetyltrimethylammonium (CTA)
molecules, their bromide ions and some water, every position folded back into
the periodic box, so the micelle is split across a face of the box. The
script

1. reads the whole file and shows the notes the reader makes,
2. reads only the micelle, with the elements from topology masses, unwrapped,
3. computes h-CTA and d-CTA in H2O and D2O in one series, and
4. shows which notes the calculation adds.

    python examples/03_micelle_contrast_series.py
"""

import warnings
from pathlib import Path

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'

# The atom masses of one CTA in topology order, as an .itp file lists them:
# N, three CH3 on the nitrogen, fifteen CH2 and the terminal CH3 of the chain.
MASS = {'N': 14.007, 'C': 12.011, 'H': 1.008}
CTA_ORDER = ['N'] + ['C', 'H', 'H', 'H'] * 3 + ['C', 'H', 'H'] * 15 + ['C', 'H', 'H', 'H']
CTA_MASSES = [MASS[element] for element in CTA_ORDER]

# 1. Everything in the file. Each note is a Diagnostic: its text, a severity
#    and the step that made it. The reader reports each once, as it reads;
#    here they are collected instead of printed.
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    sinccity.set_verbosity('quiet')
    everything = sinccity.load(DATA / 'micelle.gro')
    sinccity.set_verbosity('info')
print(f'{everything}: residues {sorted(set(everything.resnames.tolist()))}')
for note in everything.notes:
    print(f'  [{note.severity:<7} {note.origin}] {note}')
assert {note.origin for note in everything.notes} == {'read_gro'}
assert any('pass residues=' in note and note.severity == 'warning' for note in everything.notes)

# 2. The micelle only. The masses fix the elements (no guessing), and the
#    molecules are made whole and joined across the box face.
print()
gsc = GenericScatteringCalculator()
micelle = gsc.load_structure(DATA / 'micelle.gro', residues=['CTA'], masses={'CTA': CTA_MASSES})
print(f'elements from: {micelle.element_sources}')
print(f'box {micelle.box} Ang; extent after unwrapping {np.ptp(micelle.xyz, axis=0).round(1)} Ang')
print(
    f'Rg by mass {micelle.radius_of_gyration("mass").value:.2f} Ang, '
    f'by scattering length in vacuum {micelle.radius_of_gyration("sld").value:.2f} Ang'
)
assert micelle.element_sources == {'CTA': 'masses'}
assert set(micelle.resnames.tolist()) == {'CTA'} and micelle.n_points == 24 * len(CTA_MASSES)
# Folded into the box the micelle spans it; made whole it is about 50 Ang across.
assert np.all(np.ptp(micelle.xyz, axis=0) < 0.8 * np.array(micelle.box)), 'not unwrapped'

# 3. Four contrasts. Entries that share a deuteration share one pair pass.
#    The atomic volumes come from bulk elemental densities, so in a solvent
#    the absolute scale is only qualitative; the calculator says so, once per
#    entry, and records it in the profile's notes.
gsc.set_q_grid(0.01, 0.5, 25)
gsc.set_params(normalisation_volume=micelle.total_volume)
series = [
    {'label': 'h-CTA in D2O', 'params': {'sld_solvent': 6.36}},
    {'label': 'h-CTA in H2O', 'params': {'sld_solvent': -0.56}},
    {'label': 'd-CTA in H2O', 'params': {'sld_solvent': -0.56, 'deuterated': ['CTA']}},
    {'label': 'd-CTA in D2O', 'params': {'sld_solvent': 6.36, 'deuterated': ['CTA']}},
]
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # read from the notes below instead
    profiles = gsc.calculate_series(series)

columns = [0, 6, 12, 18, 24]
header = ''.join(f'{gsc.q[i]:>10.3f}' for i in columns)
print(f'\n{"I(q) in 1/cm at q =":<20}{header}{"time (s)":>10}')
for profile in profiles:
    values = ''.join(f'{profile.intensity[i]:>10.4g}' for i in columns)
    print(f'{profile.label:<20}{values}{profile.elapsed:>10.3f}')

# Hydrogenous CTA (SLD about -0.2 here) is nearly matched by H2O (-0.56), so that
# curve is weak; deuterating the chains makes them visible in H2O. The second
# entry of each deuteration state reused the first one's spectrum, hence its
# short time.
print(f'\nkernel {profiles[0].method}, bin {profiles[0].histogram_bin:.2g} Ang')
h_d2o, h_h2o, d_h2o, d_d2o = (profile.intensity[0] for profile in profiles)
assert h_h2o < 0.05 * min(h_d2o, d_h2o, d_d2o), 'h-CTA in H2O should be nearly matched'
assert d_h2o > 10 * h_h2o, 'deuteration should make the chains visible in H2O'

# 4. The notes of one profile: the structure's records first (the same
#    objects as micelle.notes), then those of the calculation.
print('\nnotes of the first profile:')
for note in profiles[0].notes:
    print(f'  [{note.severity:<7} {note.origin}] {note[:90]}{"..." if len(note) > 90 else ""}')
first = profiles[0].notes
assert first[: len(micelle.notes)] == micelle.notes, "the structure's records come first"
assert any(note.origin == 'calculate' and 'protein-like' in note for note in first)
