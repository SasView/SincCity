"""Console output: the package logger, diagnostics and encoding-safe glyphs.

Diagnostics (``sinccity.Diagnostic``, the ``notes`` of structures and
profiles) reach the console through ``emit`` alone: a ``'warning'`` as a
``UserWarning`` at the caller's line, an ``'info'`` as a log line.

Status messages that are a side effect of a call (structure loaded, q grid
set, profile computed) go through ``logging.getLogger('sinccity')`` instead of
``print``, so a caller can turn them off::

    import sinccity
    sinccity.set_verbosity('quiet')

The logger owns a stdout handler and does not propagate, so by default the
messages simply appear and do not duplicate into an application's own logging
configuration. To hand them over to that configuration instead::

    import logging
    log = logging.getLogger('sinccity')
    log.handlers.clear()
    log.propagate = True

Tables the user asks for explicitly (``describe()``) still ``print``: they
*are* the result of the call rather than a side effect of it.

The glyph constants fall back to ASCII when the active stdout cannot encode
them. A legacy Windows console (cp1252) has no ``✓`` or ``⁻``, and printing
them there raises ``UnicodeEncodeError``.
"""

import logging
import os
import sys
import warnings
from typing import Any

from ._base import Diagnostic

LOGGER_NAME = 'sinccity'

logger = logging.getLogger(LOGGER_NAME)

#: This package's directory, as frame file names spell it, ending in a separator.
_PACKAGE_DIR = os.path.join(os.path.normcase(os.path.dirname(os.path.abspath(__file__))), '')


class _StdoutHandler(logging.StreamHandler):
    """A ``StreamHandler`` that resolves ``sys.stdout`` at emit time.

    ``StreamHandler`` binds its stream at construction, which would miss every
    later redirection: the Jupyter kernel's streams, pytest's capture,
    ``contextlib.redirect_stdout``. Resolving per emit keeps the logger
    behaving like a ``print`` call.
    """

    @property
    def stream(self) -> Any:
        return sys.stdout

    @stream.setter
    def stream(self, value: Any) -> None:
        """Ignore the stream ``StreamHandler.__init__`` tries to bind."""


def _install_default_handler() -> None:
    handler = _StdoutHandler()
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False


if not logger.handlers:
    _install_default_handler()


_LEVEL_ALIASES = {'quiet': logging.WARNING}


def set_verbosity(level: int | str = 'info') -> None:
    """Set how much the calculator reports about what it is doing.

    Args:
        level: ``'quiet'`` to drop the progress messages and the ``'info'``
            notes and keep only warnings, or any standard logging level as a
            name (``'info'``, ``'debug'``, ``'warning'``, ``'error'``) or as an
            integer. ``'error'`` or above suppresses the warnings as well;
            the notes stay on the structure and the profile either way.

    Raises:
        ValueError: If *level* is not a recognized level name.
    """
    if isinstance(level, str):
        resolved = _LEVEL_ALIASES.get(level.lower(), logging.getLevelName(level.upper()))
        if not isinstance(resolved, int):
            raise ValueError(
                f"Unknown verbosity '{level}'. Use 'quiet', 'info', 'debug', "
                "'warning', 'error', or a logging level number."
            )
        level = resolved
    logger.setLevel(level)


def _caller_level() -> int:
    """The ``stacklevel`` of the first frame outside this package, seen from ``emit``.

    Counting frames instead of fixing a number keeps a warning on the user's
    line however many package frames lie between it and ``emit``.
    """
    frame, level = sys._getframe(1), 1  # level 1 is emit itself, as warnings.warn counts
    while frame is not None and os.path.normcase(frame.f_code.co_filename).startswith(_PACKAGE_DIR):
        frame, level = frame.f_back, level + 1
    return level


def emit(diagnostic: Diagnostic, context: str | None = None) -> None:
    """Report *diagnostic* once, as its severity says.

    The one place the package turns a diagnostic into output. A ``'warning'``
    is a ``UserWarning`` attributed to the first frame outside the package; an
    ``'info'`` is a log line. Both follow ``set_verbosity``.

    Args:
        diagnostic: The record.
        context: Put before the text, such as the file a reader note is about.
    """
    text = f'{context}: {diagnostic.text}' if context else diagnostic.text
    if diagnostic.severity == 'info':
        logger.info(f'  Note: {text}')
    elif logger.isEnabledFor(logging.WARNING):
        warnings.warn(text, UserWarning, stacklevel=_caller_level())


def _encodable(text: str) -> bool:
    """Return True when the current stdout can encode *text*."""
    encoding = getattr(sys.stdout, 'encoding', None) or 'utf-8'
    try:
        text.encode(encoding)
    except (LookupError, UnicodeEncodeError):
        return False
    return True


def _glyph(preferred: str, fallback: str) -> str:
    return preferred if _encodable(preferred) else fallback


#: Prefix marking a completed action in a status message.
OK = _glyph('✓', '[ok]')
#: Unit of q.
INVERSE_ANGSTROM = _glyph('Å⁻¹', '1/Ang')
#: Unit of length.
ANGSTROM = _glyph('Å', 'Ang')
#: Unit of volume.
CUBIC_ANGSTROM = _glyph('Å³', 'Ang^3')
#: Unit of scattering length density, as sasmodels parameters use it.
SLD_UNIT = _glyph('10⁻⁶ Å⁻²', '1e-6/Ang^2')
