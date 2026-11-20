"""A structure factor for a particle that is not a sphere: β(q) and the decoupling approximation.

A prolate ellipsoid of revolution (semi-axes 20, 20 and 50 Å) of SLD pixels on
a 2 Å lattice, in D2O at a 20 % volume fraction. With set_params(beta=True)
the profile also holds the orientation-averaged amplitude <F(q)> and
β(q) = <F>² / <F²>, and profile.with_structure_factor(S) applies a hard-sphere
S(q) in the decoupling approximation:

    I(q) = P(q) [1 + β(q) (S(q) - 1)]

β is checked against the continuum ellipsoid. Three curves are compared: P(q),
the plain product P(q) S(q), which treats the particle as a sphere, and the
decoupled one. β depends on the point of the particle that S(q) describes; it
is shown about the centre and about one tip.

    python examples/08_beta_structure_factor.py [--plot]
"""

import sys

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator, Structure

A, C, STEP = 20.0, 50.0, 2.0  # semi-axes and lattice step, Ang
VOLFRACTION = 0.2


def sphere_amplitude(x):
    """3 (sin x - x cos x) / x³: the normalised amplitude of a sphere at x = q R."""
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3


def ellipsoid_beta(q, a, c, n=76):
    """β of a continuum ellipsoid of revolution, averaged over cos(angle) by Gauss-Legendre."""
    mu, weight = np.polynomial.legendre.leggauss(n)
    mu, weight = (mu + 1) / 2, weight / 2  # on [0, 1]
    amplitude = sphere_amplitude(q[:, None] * np.sqrt(a**2 * (1 - mu**2) + c**2 * mu**2))
    return (amplitude @ weight) ** 2 / (amplitude**2 @ weight)


def hard_sphere(q, radius, volfraction):
    """The Percus-Yevick S(q) of hard spheres (Kinning and Thomas), as sasmodels' hardsphere."""
    phi = volfraction
    alpha = (1 + 2 * phi) ** 2 / (1 - phi) ** 4
    beta = -6 * phi * (1 + phi / 2) ** 2 / (1 - phi) ** 4
    gamma = phi * alpha / 2
    x = 2 * q * radius
    sin, cos = np.sin(x), np.cos(x)
    g = (
        alpha * (sin - x * cos) / x**2
        + beta * (2 * x * sin + (2 - x**2) * cos - 2) / x**3
        + gamma * (-(x**4) * cos + 4 * ((3 * x**2 - 6) * cos + (x**3 - 6 * x) * sin + 6)) / x**5
    )
    return 1 / (1 + 24 * phi * g / x)


# The pixels of the ellipsoid, built in memory.
axis = np.arange(-np.ceil(C / STEP), np.ceil(C / STEP) + 1) * STEP
xyz = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), -1).reshape(-1, 3)
xyz = xyz[(xyz[:, 0] ** 2 + xyz[:, 1] ** 2) / A**2 + xyz[:, 2] ** 2 / C**2 <= 1]
particle = Structure.from_arrays(xyz, np.full(len(xyz), 1e-6), np.full(len(xyz), STEP**3))

sinccity.set_verbosity('quiet')
gsc = GenericScatteringCalculator(particle)
gsc.set_q_grid(0.003, 0.3, 80)
# S(q) belongs to the number density scale / V_norm. With V_norm the particle
# volume, scale is the volume fraction S(q) is computed for.
volume = particle.total_volume
gsc.set_params(sld_solvent=6.36, scale=VOLFRACTION, normalisation_volume=volume, beta=True)
centre, tip = gsc.calculate_series(
    [
        {'label': 'about the centre'},  # beta_centre='volume', the default
        {'label': 'about a tip', 'params': {'beta_centre': (0.0, 0.0, C)}},
    ]
)
q, intensity = centre.q, centre.intensity
print(f'{particle.n_points} pixels, {volume:.0f} Ang^3; centre {centre.beta_centre_xyz} Ang')

# β against the continuum ellipsoid of the same volume as the pixels. The
# staircase surface of the pixels sets the tolerance: about 1e-2 at q_max.
k = (volume / (4 / 3 * np.pi * A**2 * C)) ** (1 / 3)
continuum = ellipsoid_beta(q, k * A, k * C)
off = np.max(np.abs(centre.beta - continuum))
print(f'beta against the continuum ellipsoid: largest difference {off:.2e} (tolerance 0.02)')
assert np.all(np.isfinite(centre.beta)) and np.all(centre.beta <= 1), 'beta lies in [0, 1]'
assert centre.beta[0] > 0.999, 'beta -> 1 at low q'
assert off < 0.02, 'the pixel ellipsoid is off the continuum one'

# A hard-sphere S(q) for the sphere of the same volume, at the same volume fraction.
radius = (3 * volume / (4 * np.pi)) ** (1 / 3)
s = hard_sphere(q, radius, VOLFRACTION)
product = intensity * s
decoupled = centre.with_structure_factor(s)
assert np.allclose(decoupled, intensity * (1 + centre.beta * (s - 1))), 'the decoupling formula'
# 0 <= beta <= 1, so the decoupled curve lies between P(q) and P(q) S(q).
assert np.all(np.abs(decoupled - intensity) <= np.abs(product - intensity) * (1 + 1e-12))

print(f'\nhard spheres of radius {radius:.1f} Ang at volume fraction {VOLFRACTION:g}')
print(
    f'{"q (1/Ang)":>10} {"S(q)":>7} {"P (1/cm)":>10} {"P S":>10} {"decoupled":>10} '
    f'{"beta":>7} {"beta tip":>8}'
)
for i in range(0, q.size, 8):
    print(
        f'{q[i]:10.4f} {s[i]:7.3f} {intensity[i]:10.4g} {product[i]:10.4g} {decoupled[i]:10.4g} '
        f'{centre.beta[i]:7.3f} {tip.beta[i]:8.3f}'
    )

# About a tip the points lie further from the centre of the amplitude, so it
# averages out at lower q: the same particle, a different β.
half = [q[np.argmax(profile.beta < 0.5)] for profile in (centre, tip)]
print(f'\nbeta falls below 0.5 at q = {half[0]:.3g} about the centre, {half[1]:.3g} about a tip')
assert half[1] < half[0] / 2, 'beta about a tip should fall at much lower q'

if '--plot' in sys.argv:
    fig = centre.plot(show=False)
    fig.data[0].name = 'P(q)'
    fig.add_scatter(x=q, y=product, mode='lines', name='P(q) S(q)')
    fig.add_scatter(x=q, y=decoupled, mode='lines', name='P(q) [1 + beta (S - 1)]')
    fig.show()
    fig = centre.plot(quantity='beta', show=False)
    fig.add_scatter(x=q, y=tip.beta, mode='lines', name=tip.label)
    fig.add_scatter(x=q, y=continuum, mode='lines', name='continuum ellipsoid')
    fig.show()
