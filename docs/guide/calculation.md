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
