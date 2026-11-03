"""Write the synthetic structure files that the scripts in ``examples/`` read.

    python tests/data/make_example_files.py

Every file is generated from fixed numbers and a fixed seed, so running this
again reproduces them byte for byte:

- ``core_shell_sphere.sld``: pixels of a core-shell sphere on a 3 Å lattice.
- ``micelle.gro``: a CTA micelle with bromide ions and some water, written
  across a face of its periodic box as GROMACS writes it.
- ``rod_hexahedra.vtk``: a 20 x 20 x 80 Å rod meshed into 5 Å hexahedra.
"""

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


# ============================================================================
# Core-shell sphere of pixels (SasView SLD file)
# ============================================================================

STEP = 3.0  # Å
CORE_RADIUS, SHELL_RADIUS = 15.0, 24.0  # Å
CORE_SLD, SHELL_SLD = 1.0e-6, 4.0e-6  # Å^-2


def core_shell_sphere() -> str:
    n = int(SHELL_RADIUS / STEP)
    axis = np.arange(-n, n + 1) * STEP
    xyz = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), -1).reshape(-1, 3)
    r = np.linalg.norm(xyz, axis=1)
    xyz, r = xyz[r <= SHELL_RADIUS], r[r <= SHELL_RADIUS]
    sld = np.where(r <= CORE_RADIUS, CORE_SLD, SHELL_SLD)
    lines = ['x y z sld']  # one header line; the pixel volume follows from the lattice
    lines += [f'{x:.1f} {y:.1f} {z:.1f} {s:.3e}' for (x, y, z), s in zip(xyz, sld)]
    return '\n'.join(lines) + '\n'


# ============================================================================
# CTA micelle (GROMACS .gro)
# ============================================================================

N_CTA = 24
BOX = 7.0  # nm, cubic
CENTRE = np.array([0.3, 3.5, 3.5])  # nm: the micelle straddles the x = 0 face
HEAD_RADIUS = 24.0  # Å, nitrogen distance from the centre
CHAIN_STEP = 1.2  # Å per CH2 along the radius (a zigzag chain projected)
N_WATER = 150
SEED = 20251103

#: Element of each CTA atom in residue order, as a topology lists them:
#: N, three methyl groups CH3, fifteen CH2 and the terminal CH3 of the chain.
CTA_ELEMENTS = ['N'] + ['C', 'H', 'H', 'H'] * 3 + ['C', 'H', 'H'] * 15 + ['C', 'H', 'H', 'H']


def _frame(u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two unit vectors perpendicular to *u* and to each other."""
    helper = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(u, helper)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(u, e1)


def _cta(u: np.ndarray) -> np.ndarray:
    """Positions (Å, about the micelle centre) of one CTA pointing along *u*."""
    e1, e2 = _frame(u)
    nitrogen = HEAD_RADIUS * u
    atoms = [nitrogen]
    for k in range(3):  # methyl groups on the outside of the nitrogen
        angle = 2 * np.pi * k / 3
        side = np.cos(angle) * e1 + np.sin(angle) * e2
        bond = 0.33 * u + 0.94 * side
        carbon = nitrogen + 1.49 * bond
        atoms.append(carbon)
        f1, f2 = _frame(bond)
        for j in range(3):  # pointing away from the nitrogen
            twist = 2 * np.pi * j / 3
            atoms.append(
                carbon + 1.09 * (0.33 * bond + 0.94 * (np.cos(twist) * f1 + np.sin(twist) * f2))
            )
    for k in range(16):  # the chain, pointing inwards
        zigzag = 0.42 if k % 2 else -0.42
        carbon = nitrogen - CHAIN_STEP * (k + 1) * u + zigzag * e1
        atoms.append(carbon)
        atoms.extend(
            [
                carbon + 0.6 * np.sign(zigzag) * e1 + 0.9 * e2,
                carbon + 0.6 * np.sign(zigzag) * e1 - 0.9 * e2,
            ]
        )
        if k == 15:
            atoms.append(carbon - 1.09 * u)
    return np.array(atoms)


def _cta_names() -> list[str]:
    names, carbon, hydrogen = ['N1'], 0, 0
    for element in CTA_ELEMENTS[1:]:
        if element == 'C':
            carbon += 1
            names.append(f'C{carbon}')
        else:
            hydrogen += 1
            names.append(f'H{hydrogen}')
    return names


def _directions(n: int) -> np.ndarray:
    """*n* nearly even unit vectors (a Fibonacci sphere)."""
    k = np.arange(n) + 0.5
    polar = np.arccos(1 - 2 * k / n)
    azimuth = np.pi * (1 + 5**0.5) * k
    return np.column_stack(
        [np.sin(polar) * np.cos(azimuth), np.sin(polar) * np.sin(azimuth), np.cos(polar)]
    )


def micelle() -> str:
    rng = np.random.default_rng(SEED)
    residues = []  # (resname, [(name, xyz in nm, about the centre)])
    names = _cta_names()
    for u in _directions(N_CTA):
        residues.append(('CTA', list(zip(names, _cta(u) / 10))))
    for u in _directions(N_CTA) @ _random_rotation(rng).T:  # one bromide per head group
        residues.append(('BR', [('BR', 2.95 * u)]))
    placed = 0
    while placed < N_WATER:
        oxygen = rng.uniform(-BOX / 2, BOX / 2, 3)
        if np.linalg.norm(oxygen) < 3.3:  # nm: outside the micelle and its ions
            continue
        e1, e2 = _frame(oxygen / np.linalg.norm(oxygen))
        hydrogens = [oxygen + 0.09572 * (np.cos(a) * e1 + np.sin(a) * e2) for a in (0.0, 1.824)]
        residues.append(('SOL', [('OW', oxygen), ('HW1', hydrogens[0]), ('HW2', hydrogens[1])]))
        placed += 1

    lines, number = ['CTA micelle, 24 CTA + 24 BR + 150 SOL (synthetic, for the examples)'], 0
    body = []
    for resid, (resname, atoms) in enumerate(residues, start=1):
        for name, xyz in atoms:
            number += 1
            x, y, z = (CENTRE + xyz) % BOX  # every position back in the box, as GROMACS writes it
            body.append(f'{resid:>5d}{resname:<5}{name:>5}{number:>5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    lines.append(f'{number:5d}')
    lines += body
    lines.append(f'{BOX:10.5f}{BOX:10.5f}{BOX:10.5f}')
    return '\n'.join(lines) + '\n'


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.diag(r))


# ============================================================================
# Rod of hexahedra (legacy ASCII VTK)
# ============================================================================

ROD = (20.0, 20.0, 80.0)  # Å
CELL = 5.0  # Å
ROD_SLD = 3.0e-6  # Å^-2


def rod_hexahedra() -> str:
    counts = [int(round(side / CELL)) for side in ROD]
    axes = [np.arange(n + 1) * CELL - side / 2 for n, side in zip(counts, ROD)]
    points = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    nx, ny, nz = (n + 1 for n in counts)

    def index(i, j, k):
        return (i * ny + j) * nz + k

    corners = (
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    )
    cells = [
        [index(i + a, j + b, k + c) for a, b, c in corners]
        for i in range(counts[0])
        for j in range(counts[1])
        for k in range(counts[2])
    ]
    lines = [
        '# vtk DataFile Version 3.0',
        'Rod 20 x 20 x 80 Ang in 5 Ang hexahedra (synthetic, for the examples)',
        'ASCII',
        'DATASET UNSTRUCTURED_GRID',
        f'POINTS {len(points)} float',
        *(f'{x:g} {y:g} {z:g}' for x, y, z in points),
        f'CELLS {len(cells)} {9 * len(cells)}',
        *(' '.join(map(str, [8, *cell])) for cell in cells),
        f'CELL_TYPES {len(cells)}',
        *(['12'] * len(cells)),
        f'CELL_DATA {len(cells)}',
        'SCALARS sld float 1',
        'LOOKUP_TABLE default',
        *([f'{ROD_SLD:g}'] * len(cells)),
    ]
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    for name, text in (
        ('core_shell_sphere.sld', core_shell_sphere()),
        ('micelle.gro', micelle()),
        ('rod_hexahedra.vtk', rod_hexahedra()),
    ):
        (HERE / name).write_text(text, encoding='ascii', newline='\n')
        print(f'wrote {name}: {text.count(chr(10))} lines')
