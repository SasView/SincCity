"""Reproducible settings records, replayed in other processes.

A settings record names the file (relative to a root, so it survives a move
of the data), its SHA-256, the reader options, the q grid, every parameter
and the library versions. It is plain JSON, and so is a series of contrasts,
so a task can cross a process boundary or be stored next to the results.
This script saves a record, replays it in two worker processes that split
a series between them, and checks the curves against the in-process ones.

    python examples/07_settings_and_workers.py
"""

import json
import multiprocessing
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sinccity
from sinccity import GenericScatteringCalculator

DATA = Path(__file__).resolve().parents[1] / 'tests' / 'data'
OUTPUT = Path(__file__).resolve().parent / 'output'

SERIES = [
    {'label': 'h-CTA in D2O', 'params': {'sld_solvent': 6.36}},
    {'label': 'h-CTA in H2O', 'params': {'sld_solvent': -0.56}},
    {'label': 'd-CTA in H2O', 'params': {'sld_solvent': -0.56, 'deuterated': ['CTA']}},
    {'label': 'd-CTA in D2O', 'params': {'sld_solvent': 6.36, 'deuterated': ['CTA']}},
]


def run(task: str) -> dict[str, list[float]]:
    """What a worker runs: a JSON task of settings, series and data root."""
    sinccity.set_verbosity('quiet')  # the logger is per process
    task = json.loads(task)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # the notes travel back with the profiles if needed
        calculator = GenericScatteringCalculator.from_settings(task['settings'], root=task['root'])
        profiles = calculator.calculate_series(task['series'])
    return {profile.label: profile.intensity.tolist() for profile in profiles}


def main() -> None:
    sinccity.set_verbosity('quiet')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # elements guessed from names: fine for these atoms
        gsc = GenericScatteringCalculator()
        gsc.load_structure(DATA / 'micelle.gro', residues=['CTA'])
        gsc.set_q_grid(0.01, 0.5, 40)
        gsc.set_params(normalisation_volume=gsc.structure.total_volume, size_pd=0.1)
        expected = {p.label: p.intensity for p in gsc.calculate_series(SERIES)}

    # The path is stored relative to the stated root, with '/' separators.
    settings = gsc.get_settings(root=DATA)
    OUTPUT.mkdir(exist_ok=True)
    record = OUTPUT / 'micelle-settings.json'
    record.write_text(json.dumps(settings, indent=2))
    print(f'saved {record}')
    print(f'  source {settings["source"]!r}, sha256 {settings["source_sha256"][:12]}...')
    print(f'  reader options {settings["reader_options"]}')
    print(f'  versions {settings["versions"]}')

    # Two workers, each with half the series; 'spawn' starts clean interpreters,
    # as on Windows and macOS. The root is where the data lives for the workers.
    saved = json.loads(record.read_text())
    tasks = [
        json.dumps({'settings': saved, 'series': half, 'root': str(DATA)})
        for half in (SERIES[:2], SERIES[2:])
    ]
    context = multiprocessing.get_context('spawn')
    with ProcessPoolExecutor(max_workers=2, mp_context=context) as pool:
        results = {}
        for part in pool.map(run, tasks):
            results.update(part)

    print(f'\n{"curve":<14} {"I(q_min) (1/cm)":>16} {"largest difference from in-process":>36}')
    for label, intensity in results.items():
        difference = np.max(np.abs(np.array(intensity) / expected[label] - 1))
        print(f'{label:<14} {intensity[0]:16.5g} {difference:36.1e}')


if __name__ == '__main__':  # required for 'spawn': workers import this file
    main()
