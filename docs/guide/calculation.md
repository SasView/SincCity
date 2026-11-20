# Calculation

## The Debye sum

`calculate()` sums over all points:

$$
\begin{aligned}
w_j &= (\rho_j - \rho_s)\, V_j \\
D(q) &= \sum_j \sum_k w_j w_k \, \frac{\sin(q r_{jk})}{q r_{jk}} \\
I(q) &= \mathrm{scale} \cdot 10^8 \, \frac{D(q)}{V_\mathrm{norm}} + \mathrm{background}
\end{aligned}
$$

with $\rho_j$ the SLD of point $j$ in Å⁻², $\rho_s$ = `sld_solvent` × 10⁻⁶
Å⁻², $V_j$ its volume in Å³ and $r_{jk}$ in Å. The factor $10^8$ turns Å⁻¹
into cm⁻¹, the unit of $I$.

`sld_solvent` is in 10⁻⁶ Å⁻², like every sasmodels SLD parameter. SasView's
calculator takes Å⁻² instead; a value that looks like Å⁻² triggers a warning.
The structure's own SLDs stay in Å⁻², as the files store them.

## Normalisation

$V_\mathrm{norm}$ is the volume of the points that scatter, those whose SLD differs
from the solvent's. For a particle of pixels or cells in a box of solvent, the
curve at `scale=1` then equals sasmodels' form factor at `scale=1`, so `scale`
is a volume fraction.

When the solvent exactly matches a population of points, those points stop
counting and the curve jumps by the volume ratio; for contrast-variation
series set `set_params(normalisation_volume=...)` explicitly (`'auto'`
restores the default). SasView calls the override `total_volume` and, when
some points have zero contrast, also multiplies by the ratio of all points'
volume to the scattering volume. This calculator does not.

## Size polydispersity

`set_params(size_pd=0.1)` sums blown-up and shrunk copies of the whole
structure, as sasmodels does for `radius_pd`: `size_pd_n` size factors (35)
evenly over ±`size_pd_nsigma` (3) standard deviations, factors at or below 0
dropped, each weighted by the Gaussian's height there. A factor $f$ stretches
every distance by $f$ and every amplitude by $f^3$, so

$$
I(q) = \mathrm{scale} \cdot 10^8 \, \frac{\sum_f n_f f^6 D(f q)}{V_\mathrm{norm}}
+ \mathrm{background}, \qquad V_\mathrm{norm} = \sum_f n_f f^3 V
$$

with $V$ the volume of the points that scatter and the shares $n_f$ adding up
to 1, so `scale` stays a volume fraction. The lattice, histogram and average
kernels evaluate their one distance spectrum at $f q$ for every size, which costs only the sums over q; `'pairs'` repeats the
whole pair sum per size and warns. The automatic histogram bin and the
cell-size warning use the largest factor, `1 + size_pd_nsigma * size_pd`.

## Deuteration

`set_params(deuterated=['CTA'])` gives the hydrogens of these residue names
the SLD of D at `calculate()`, as
[`Structure.deuterated`](structures.md#deuteration) does.

## Kernels

`set_params(kernel=...)` picks how D(q) is summed; the profile's `method`
records it.

| `kernel` | How | Exact? | Time |
|---|---|---|---|
| `'lattice'` | FFT autocorrelation of the pixel grid | yes | well under a second for 10⁵ pixels |
| `'histogram'` | pairs counted by distance, then one sum per q over the occupied distances | to about 1e-4 | about 2 s for 11 000 atoms, N² |
| `'average'` | ⟨A⟩²: each point's distance from the bounding-box centre | only for spherically symmetric structures | linear in N |
| `'pairs'` | every pair at every q | to the last digit | 2–3 s for 4 000 points, 25 s for 8 000 |

`'auto'` (the default) is `'lattice'` for pixels on a lattice and
`'histogram'` for everything else.

The histogram bin is `'auto'`:
$\min\left(0.01~\text{Å},\ 10^{-3} / (q_\mathrm{max} f_\mathrm{max})\right)$,
with $f_\mathrm{max}$ the largest size factor (1 without polydispersity),
which keeps the rounding phase error below $5 \cdot 10^{-4}$ rad; `histogram_bin=` sets it, and
`profile.histogram_bin` records it. A cached spectrum with a finer bin is
reused.

`'average'` reproduces SasView's `is_avg=True` and suits round aggregates
above about 10⁵ atoms. `calculate(progress=...)` can report progress or stop
the pair passes of `'histogram'` and `'pairs'`.

The amplitudes of `'average'` are taken about the centre of the bounding box,
the origin SasView's `transform_center` uses, so the two agree. For a
lopsided aggregate the centre of scattering length would be the more natural
origin, and the result depends on the choice; the profile's note records the
centre used, in Å.

`'lattice'` zero-pads the grid to $2n - 1$ cells per axis and refuses grids
above 8 million padded cells (`sinccity.debye.LATTICE_MAX_CELLS`); `'auto'`
then falls back to `'histogram'`. The cap counts cells, not bytes: at the cap
the transforms and their temporaries peak near 400 MB.

## β(q) and structure factors

`set_params(beta=True)` adds the orientation-averaged amplitude about a centre
$c$ and the ratio β(q) to the profile:

$$
A(q) = \sum_j w_j \operatorname{sinc}(q |x_j - c|), \qquad
\beta(q) = \frac{A(q)^2}{D(q)}
$$

$D - A^2$ is the variance of the amplitude over orientations, so
$0 \le \beta \le 1$. β = 1 at every q only for a spherically symmetric
particle about its centre. As q → 0, β → 1 while the net contrast $\sum w$ is
not zero; at zero net contrast β falls as q² instead, and close to it β drops
from 1 within a q range that shrinks with $|\sum w|$.

`profile.amplitude` (Å, signed), `profile.beta` and `profile.beta_centre_xyz`
(Å) hold them. `I(q)` itself does not change. The amplitude costs time linear
in the points, so it adds little to the pair pass, and the spectrum cache
works as before. It cannot be combined with `kernel='average'`, whose P(q) is
itself ⟨A⟩²: a series entry that sets `kernel='average'` on a calculator with
`beta=True` must also set `beta=False`.

### Applying S(q)

For particles whose orientations are independent of their positions, the
decoupling approximation gives

$$
I(q) = \mathrm{scale} \cdot 10^8 \, \frac{D(q) + A(q)^2 \left(S(q) - 1\right)}{V_\mathrm{norm}}
+ \mathrm{background}
= P(q) \left[1 + \beta(q) \left(S(q) - 1\right)\right] + \mathrm{background}
$$

`profile.with_structure_factor(s)` returns it, for *s* an array with one value
per q or a callable of q. It uses the first form, which never divides by
D(q), so it is finite where β is nan.

```python
gsc.set_params(sld_solvent=6.36, beta=True,
               scale=0.2, normalisation_volume=particle_volume)
profile = gsc.calculate()
intensity = profile.with_structure_factor(hard_sphere(profile.q, radius, 0.2))
```

[Example 08](https://github.com/SasView/SincCity/tree/master/examples/08_beta_structure_factor.py)
writes out a Percus–Yevick `hard_sphere` and checks β of a pixel ellipsoid
against the continuum one.

!!! warning "What the approximation assumes"
    1. **The structure is one whole particle.** D(q) does not need that; A(q)
       about one centre does. Several aggregates, a micelle with free
       surfactant, or an aggregate split across a periodic face give a β that
       looks plausible and means nothing. The profile's info note gives the
       largest and the rms distance of the points from the centre, which shows
       such a case as a size that does not fit.
    2. **Orientations are independent of positions.** Dense or strongly
       elongated particles break this.
    3. **Sizes are independent of positions.** With `size_pd > 0` one S(q)
       serves every size; for excluded-volume interactions it must then be an
       effective structure factor for the distribution.
    4. **S(q) belongs to the number density n = scale / V_norm.** Its volume
       fraction is n times the volume the particles exclude. `scale` is that
       fraction only when V_norm is that volume. It is not when a component is
       matched (V_norm is then the rest, and a note says so) or for atoms,
       whose elemental volumes are not excluded volumes. Set
       `normalisation_volume` to the particle volume, always in a contrast
       series that S(q) is applied to.

### The centre

S(q) describes where the particle centres are, so the centre is the point of
each particle that S(q) is about. D(q) does not depend on it; A(q) and β do.

| `beta_centre` | Centre | Moves within a contrast series |
|---|---|---|
| `'volume'` (default) | volume-weighted centroid of all points | no |
| `'active'` | volume-weighted centroid of the points that scatter | when a population matches exactly |
| `'contrast'` | centroid of contrast × volume, SasView's in vacuum | yes |
| `'box'` | middle of the bounding box, as `kernel='average'` uses | no |
| `(x, y, z)` | as given, Å | no |

β changes through a contrast series because the weights change; the default
keeps the point it is taken about fixed. `'active'` suits a pixel file with a
box of solvent around an off-centre particle (when the particle sits in the
middle of its box, `'volume'` and `'box'` find it too). `'contrast'` is an
error at zero net contrast, and when the centroid falls outside the points
that scatter, which happens close to it. An explicit centre is in the
coordinates of `structure.xyz`, after any unwrapping by the reader; a warning
says when it lies outside the points that scatter.

### Sizes

A size factor f contributes f³ A(fq) to the amplitude, which is squared after
the sum over sizes:

$$
\beta(q) = \frac{\left(\sum_f n_f f^3 A(f q)\right)^2}{\sum_f n_f f^6 D(f q)}
$$

This is sasmodels' β, F1² / F2 over the distribution. β ≤ 1 still holds
(Cauchy–Schwarz), and with a spread β → (Σ n f³)² / Σ n f⁶ < 1 as q → 0. Each
size is a copy blown up about the centre.

### When β is nan

β is a ratio, so it is unreliable wherever D(q) is small against its error:
at the minima of P(q), and at low q near zero net contrast. β is nan where
D(q) is below 10 times an estimate of that error, and a warning note gives
how many points and their q range. For `'histogram'` the estimate is the
largest phase error of a pair, $q f_\mathrm{max} \cdot \mathrm{bin} / 2$,
times $\sum_j w_j^2$, the pair weights summed in quadrature; for the exact
kernels it is a round-off allowance of $10^{-10} (\sum_j |w_j|)^2$. The
quadrature sum assumes the rounding errors of different pairs are
independent, which holds for atoms from a simulation; a structure with many
pairs at exactly one distance, such as a lattice forced through
`'histogram'`, can have a larger error. A β above 1 by more than 10⁻³ also
gives a warning.

### Output

`profile.to_data1d(quantity='beta')`, `profile.save_csv(path,
quantity='beta')` and `profile.plot(quantity='beta')` give β(q) without its
nan points, and say how many were left out. A β plot has a linear y axis and
takes no measured data. `profile.describe()` gives the range of β and the
centre.

SasView's Generic Scattering Calculator has a β(Q) option too. Its amplitude
uses the scattering lengths in vacuum about the centroid of scattering length,
for atoms only, and divides β by its first value. Here the amplitude uses the
same contrast as P(q) and needs no rescaling; `beta_centre='contrast'` with
`sld_solvent=0` gives SasView's shape.

## Things to know

!!! warning "PDB structures in a solvent are only qualitative"
    The natural atomic volumes use bulk elemental densities, liquid H₂, N₂
    and O₂ for H, N and O. A typical residue gets about 310 Å³ instead of
    about 133 Å³, so in D₂O a protein comes out about 8× too intense at
    q → 0. There is no hydration shell and no H/D exchange either.
    `calculate()` warns with the actual volume ratio. In vacuum
    (`sld_solvent=0`) the result is right. For atoms, `scale` is not a
    volume fraction.

- **Points stand in for cells.** Pixels and VTK cells are treated as points,
  which degrades as q × cell size approaches 1. `calculate()` warns above
  that, with the cell size taken as the lattice step or the cell diameter
  divided by √3.
- **The curve is theory.** `profile.to_data1d()` returns a `Data1D` with no
  dI and no dQ, with `qmin`, `qmax` and an empty `mask` set. A fit of it
  should state its weighting and resolution explicitly; pass
  `to_data1d(dI=...)` to attach uncertainties.
- **Notes.** What a reader dropped or assumed and what a calculation warns
  about are `Diagnostic` records in `structure.notes` and `profile.notes`:
  each is its text as a string, with a `severity` (`'warning'` or `'info'`)
  and an `origin` (`'read_pdb'`, `'load'`, `'calculate'`, ...). Each is
  reported once, where it is made: a reader's notes when the file is read,
  not again when the structure is installed. Warnings are `UserWarning`s
  pointing at your line; info notes are log lines. `profile.notes` holds the
  structure's records followed by the calculation's.
- **Console output.** Progress lines and info notes go through the `sinccity`
  logger; `sinccity.set_verbosity('quiet')` silences them and
  `set_verbosity('error')` the warnings too.
