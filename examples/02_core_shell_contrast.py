"""Contrast variation of a core-shell particle given as SLD pixels.

core_shell_sphere.sld holds the pixels of a sphere on a 3 Å lattice: a core
of radius 15 Å (SLD 1e-6 Å⁻²) in a shell to 24 Å (4e-6 Å⁻²). Pixels on a
lattice use the exact FFT kernel. One series computes the particle in four
solvents; matching the shell leaves the core alone, which is checked against
the sphere form factor, and every curve is checked against the core-shell
form factor.

    python examples/02_core_shell_contrast.py
"""

from pathlib import Path

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'


def sphere_amplitude(q, radius):
    """3 (sin x - x cos x) / x³ at x = q R: the normalised amplitude of a sphere."""
    x = q * radius
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3


def core_shell(q, core, outer, sld_core, sld_shell, sld_solvent, v_norm):
    """I(q) in 1/cm of a core-shell sphere, normalised by *v_norm* as the calculator does."""
    v_core, v_outer = 4 / 3 * np.pi * core**3, 4 / 3 * np.pi * outer**3
    amplitude = (sld_core - sld_shell) * v_core * sphere_amplitude(q, core)
    amplitude = amplitude + (sld_shell - sld_solvent) * v_outer * sphere_amplitude(q, outer)
    return 1e8 * amplitude**2 / v_norm


sinccity.set_verbosity('quiet')  # keep the output to the tables below
gsc = GenericScatteringCalculator(DATA / 'core_shell_sphere.sld')
structure = gsc.structure
print(structure)
print(f'lattice step {structure.lattice.step} Ang, pixel volume {structure.volume[0]:g} Ang^3')

# The pixel sphere is not quite a sphere: compare with spheres of the same
# volume as the core pixels and as all pixels.
core = structure.sld < 2e-6
r_core = (3 * structure.volume[core].sum() / (4 * np.pi)) ** (1 / 3)
r_outer = (3 * structure.total_volume / (4 * np.pi)) ** (1 / 3)
print(f'equivalent radii: core {r_core:.2f} Ang (nominal 15), outer {r_outer:.2f} Ang (nominal 24)')

# For a contrast series V_norm should not jump when a population is matched,
# so fix it at the particle volume. sld_solvent is in 1e-6 Ang^-2, as in sasmodels.
gsc.set_q_grid(0.005, 0.3, 40)
gsc.set_params(normalisation_volume=structure.total_volume)
solvents = {'H2O': -0.56, 'core matched': 1.0, 'shell matched': 4.0, 'D2O': 6.36}
profiles = gsc.calculate_series(
    [{'label': name, 'params': {'sld_solvent': sld}} for name, sld in solvents.items()]
)

print(
    f'\n{"solvent":<14} {"kernel":<8} {"I(q_min) (1/cm)":>16} {"vs core-shell model, q < 0.15":>31}'
)
for profile, (name, sld) in zip(profiles, solvents.items()):
    model = core_shell(profile.q, r_core, r_outer, 1e-6, 4e-6, sld * 1e-6, structure.total_volume)
    low = profile.q < 0.15
    ratio = profile.intensity[low] / model[low]
    print(
        f'{name:<14} {profile.method:<8} {profile.intensity[0]:16.4g} '
        f'{"ratio " + f"{ratio.min():.3f} to {ratio.max():.3f}":>31}'
    )

# With the shell matched only the core scatters: a sphere of contrast -3e-6.
matched = profiles[2]
sphere = 1e8 * (3e-6 * structure.volume[core].sum()) ** 2 * sphere_amplitude(matched.q, r_core) ** 2
sphere /= structure.total_volume
print('\nshell matched, against a sphere of the core volume:')
print(f'{"q (1/Ang)":>10} {"I(q)":>12} {"sphere":>12} {"ratio":>7}')
for i in range(0, matched.q.size, 6):
    print(
        f'{matched.q[i]:10.4f} {matched.intensity[i]:12.5g} {sphere[i]:12.5g} '
        f'{matched.intensity[i] / sphere[i]:7.3f}'
    )
print(f'\n{matched.n_zero_contrast} of {structure.n_points} pixels have zero contrast')
