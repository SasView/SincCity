# Examples

Each script is self-contained and runs from the repository root:

```bash
pip install -e ".[plot,sasdata]"   # plotly and sasdata are optional
python examples/01_quickstart_pdb.py
```

The structure files are in [`tests/data`](../tests/data). Files the scripts
write (CSV curves, a settings record) go to `examples/output/`, which git
ignores.

| Script | Task | Data | Shows |
|---|---|---|---|
| [`01_quickstart_pdb.py`](01_quickstart_pdb.py) | I(q) of C60 in vacuum | `c60.pdb` | `load_structure` with `elements=`, `describe`, the histogram and pair kernels, a check against a spherical shell, `save_csv`; `--plot` draws the points and the curve |
| [`02_core_shell_contrast.py`](02_core_shell_contrast.py) | Contrast variation of a core-shell particle | `core_shell_sphere.sld` | Pixels on a lattice (exact FFT kernel), `calculate_series` over four solvents with a fixed `normalisation_volume`, contrast matching, checks against the analytic core-shell and sphere form factors |
| [`03_micelle_contrast_series.py`](03_micelle_contrast_series.py) | A CTA micelle from an MD frame | `micelle.gro` | Reader notes as `Diagnostic` records, `residues=` and `masses=`, unwrapping across the periodic box, deuteration, a contrast series sharing pair passes, the atomic-volume caveat |
| [`04_vtk_mesh.py`](04_vtk_mesh.py) | A finite-element mesh | `rod_hexahedra.vtk`, `five_tetrahedra_cube.vtk` | One point per cell, the cell-size warning, a check against the orientation-averaged box form factor, SLD given per point |
| [`05_structures_in_memory.py`](05_structures_in_memory.py) | Structures without a file | none | `Structure.from_arrays` (a Janus dumbbell) and the four kernels compared, `Structure.from_atoms` (a benzene cluster) with H and D |
| [`06_polydispersity_and_export.py`](06_polydispersity_and_export.py) | A size distribution, and output | `core_shell_sphere.sld` | `size_pd` filling in the form-factor minimum, CSV, a sasdata `Data1D` with uncertainties; `--plot` overlays the curves |
| [`07_settings_and_workers.py`](07_settings_and_workers.py) | Reproducible tasks in other processes | `micelle.gro` | `get_settings(root=)` saved as JSON, `from_settings(root=)` in two `spawn` worker processes that split a series |

The three synthetic files (`core_shell_sphere.sld`, `micelle.gro`,
`rod_hexahedra.vtk`) are written by
[`tests/data/make_example_files.py`](../tests/data/make_example_files.py).
`tests/test_examples.py` runs every script, so they stay in step with the
package.
