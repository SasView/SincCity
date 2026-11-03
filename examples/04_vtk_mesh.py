"""A finite-element mesh from a legacy VTK file: one point per cell.

rod_hexahedra.vtk is a 20 x 20 x 80 Å rod meshed into 5 Å hexahedra with an
SLD per cell. Each cell becomes a point at its volume centroid, with the
cell's volume. The curve is compared with the orientation-averaged form
factor of a rectangular box, which it follows until q times the cell size
nears 1; the calculator warns past that. five_tetrahedra_cube.vtk shows an
SLD given per point, which the reader averages to the cells.

    python examples/04_vtk_mesh.py
"""

import warnings
from pathlib import Path

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'


def box_intensity(q, a, b, c, contrast, nodes=120):
    """I(q) in 1/cm of a rectangular box a x b x c (Å), averaged over orientations."""
    x, w = np.polynomial.legendre.leggauss(nodes)
    mu = (x + 1) / 2  # cos(theta) in [0, 1]
    phi = (x + 1) * np.pi / 4  # phi in [0, pi/2]
    weight = w / 2  # each set of weights sums to 1: an average, not an integral
    sin_theta = np.sqrt(1 - mu**2)
    qx = q[:, None, None] * sin_theta[None, :, None] * np.cos(phi)[None, None, :]
    qy = q[:, None, None] * sin_theta[None, :, None] * np.sin(phi)[None, None, :]
    qz = q[:, None, None] * mu[None, :, None]
    shape = np.sinc(qx * a / 2 / np.pi) * np.sinc(qy * b / 2 / np.pi) * np.sinc(qz * c / 2 / np.pi)
    average = np.einsum('qij,i,j->q', shape**2, weight, weight)
    return 1e8 * contrast**2 * a * b * c * average


sinccity.set_verbosity('quiet')
gsc = GenericScatteringCalculator(DATA / 'rod_hexahedra.vtk')
rod = gsc.structure
print(rod)
print(
    f'{rod.n_points} cells, total volume {rod.total_volume:g} Ang^3, '
    f'largest cell diameter {rod.cell_diameter.max():.3f} Ang'
)
rg_box = np.sqrt((20**2 + 20**2 + 80**2) / 12)
print(
    f'Rg of the cell centroids {rod.radius_of_gyration().value:.2f} Ang, of the box {rg_box:.2f} Ang'
)

# In vacuum the contrast is the SLD itself, 3e-6 1/Ang^2. The grid goes past
# q = 1 / cell size, where one point per cell stops being a good model.
gsc.set_q_grid(0.005, 0.5, 30)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    profile = gsc.calculate()
print(f'\nwarned: {caught[0].message}' if caught else '\nno warning')

model = box_intensity(profile.q, 20.0, 20.0, 80.0, 3e-6)
print(f'\n{"q (1/Ang)":>10} {"I(q) (1/cm)":>13} {"box":>13} {"ratio":>7}')
for i in range(0, profile.q.size, 3):
    print(
        f'{profile.q[i]:10.4f} {profile.intensity[i]:13.5g} {model[i]:13.5g} '
        f'{profile.intensity[i] / model[i]:7.3f}'
    )

# An SLD per point is averaged to each cell, and a magnetic vector field is
# ignored; both are info notes on the structure.
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cube = sinccity.load(DATA / 'five_tetrahedra_cube.vtk')
print(f'\n{cube}: volumes {cube.volume.round(4).tolist()} Ang^3 (sum {cube.total_volume:g})')
for note in cube.notes:
    print(f'  [{note.severity}] {note}')
