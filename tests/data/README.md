# Test structure files

Used by the tests in `tests/` and the scripts in `examples/`.

| File | Content | Origin |
|---|---|---|
| `c60.pdb` | Buckminsterfullerene, 60 carbon atoms, one pentagon of `CONECT` bonds | SasView test data |
| `five_tetrahedra_cube.vtk` | The cube [-1, 1]³ as five tetrahedra; nuclear SLD 1 as point data, a magnetic vector per cell | SasView test data |
| `core_shell_sphere.sld` | 2109 pixels on a 3 Å lattice: a core of radius 15 Å at SLD 1e-6 Å⁻² in a shell to 24 Å at 4e-6 Å⁻² | `make_example_files.py` |
| `micelle.gro` | 24 CTA molecules (62 atoms each), 24 Br⁻ and 150 waters in a 7 nm box, folded into the box so the micelle straddles the x = 0 face | `make_example_files.py` |
| `rod_hexahedra.vtk` | A 20 × 20 × 80 Å rod as 256 hexahedra of 5 Å, SLD 3e-6 Å⁻² per cell | `make_example_files.py` |

The two SasView files are copied unchanged from SasView
(<https://github.com/SasView/sasview>, `test/sascalculator/data/`),
distributed under the BSD 3-Clause licence, Copyright (c) 2009-2026, SasView Developers.

The other three are synthetic. `python tests/data/make_example_files.py`
rewrites them byte for byte from fixed numbers and a fixed seed.
