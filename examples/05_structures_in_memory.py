"""Structures built in memory: pixels from arrays and atoms from another program.

No file is needed. Part 1 builds a dumbbell of two pixel spheres with
``Structure.from_arrays`` and compares the four kernels: the exact lattice
FFT, the distance histogram, the exact pair sum and the ⟨A⟩² average, which
is exact only for spherically symmetric shapes. Part 2 builds a cluster of
benzene molecules with ``Structure.from_atoms``, as coordinates from an MD
analysis package would arrive, and compares it with its deuterated copy.

    python examples/05_structures_in_memory.py
"""

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator, Structure
from sinccity.elements import constants

sinccity.set_verbosity('quiet')
rng = np.random.default_rng(1)

# ----------------------------------------------------------------------------
# 1. Pixels: two touching spheres on a 2 Ang lattice
# ----------------------------------------------------------------------------

step, radius = 2.0, 12.0
axis = np.arange(-17, 18) * step
grid = np.stack(np.meshgrid(axis, axis[:13], axis[:13], indexing='ij'), -1).reshape(-1, 3)
grid[:, 1:] -= axis[6]  # centre the short axes on 0
left = np.linalg.norm(grid - [-radius, 0, 0], axis=1) <= radius
right = np.linalg.norm(grid - [radius, 0, 0], axis=1) <= radius
inside = left | right
xyz = grid[inside]
sld = np.where(left[inside], 2e-6, 5e-6)  # a Janus dumbbell, SLD in Ang^-2
dumbbell = Structure.from_arrays(xyz, sld, np.full(len(xyz), step**3))
print(f'{dumbbell}, lattice step {dumbbell.lattice.step} Ang')

gsc = GenericScatteringCalculator(dumbbell)
gsc.set_q_grid(0.01, 0.4, 30)
gsc.set_params(sld_solvent=1.0)
results = {}
for kernel in ('lattice', 'histogram', 'pairs', 'average'):
    gsc.set_params(kernel=kernel)
    results[kernel] = gsc.calculate()

# Differences on the scale of the curve, I(0) of the exact lattice kernel.
reference = results['lattice'].intensity
print(f'\n{"kernel":<10} {"time (s)":>9} {"largest difference from lattice, of I(0)":>42}')
for kernel, profile in results.items():
    difference = np.max(np.abs(profile.intensity - reference)) / reference[0]
    print(f'{kernel:<10} {profile.elapsed:9.3f} {difference:42.1e}')
print('the <A>^2 kernel is far off: a dumbbell is not spherically symmetric')

# ----------------------------------------------------------------------------
# 2. Atoms: 27 benzene molecules, hydrogenated and deuterated
# ----------------------------------------------------------------------------

angles = np.arange(6) * np.pi / 3
ring = np.column_stack([np.cos(angles), np.sin(angles), np.zeros(6)])
benzene = np.vstack([1.39 * ring, 2.48 * ring])  # six C, then six H
symbols, positions, resids, names = [], [], [], []
centres = np.stack(np.meshgrid(*[np.arange(3) * 7.0] * 3, indexing='ij'), -1).reshape(-1, 3)
for number, centre in enumerate(centres, start=1):
    rotation, upper = np.linalg.qr(rng.normal(size=(3, 3)))
    positions.append(benzene @ (rotation * np.sign(np.diag(upper))).T + centre)
    symbols += ['C'] * 6 + ['H'] * 6
    resids += [number] * 12
    names += [f'C{k}' for k in range(1, 7)] + [f'H{k}' for k in range(1, 7)]
cluster = Structure.from_atoms(
    np.vstack(positions), symbols, resids=resids, resnames=['BNZ'] * len(symbols), names=names
)
print(f'\n{cluster}, element sources {cluster.element_sources}')

gsc = GenericScatteringCalculator(cluster)
gsc.set_q_grid(0.005, 0.5, 20)
h = gsc.calculate(label='h-benzene')
gsc.set_params(deuterated=['BNZ'])  # or calculate on cluster.deuterated(['BNZ'])
d = gsc.calculate(label='d-benzene')

# In vacuum I(0) is proportional to the squared total scattering length.
b = {element: constants(element)[0] * constants(element)[1] for element in ('C', 'H', 'D')}
expected = ((b['C'] + b['H']) / (b['C'] + b['D'])) ** 2
print(f'I(0) h / d: computed {h.intensity[0] / d.intensity[0]:.4f}, expected {expected:.4f}')
print(
    f'Rg by mass {cluster.radius_of_gyration("mass").value:.2f} Ang; by scattering length: '
    f'h {cluster.radius_of_gyration().value:.2f} Ang, '
    f'd {cluster.deuterated("BNZ").radius_of_gyration().value:.2f} Ang'
)
