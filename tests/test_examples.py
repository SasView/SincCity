"""Every script in examples/ runs to the end."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'
SCRIPTS = sorted(EXAMPLES.glob('[0-9][0-9]_*.py'))


def test_there_are_examples():
    assert len(SCRIPTS) >= 7


@pytest.mark.parametrize('script', SCRIPTS, ids=[script.stem for script in SCRIPTS])
def test_example_runs(script):
    env = dict(os.environ, PYTHONIOENCODING='utf-8')  # the glyphs, on any console
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        encoding='utf-8',
        env=env,
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-2000:]
