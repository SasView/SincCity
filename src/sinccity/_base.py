"""Definitions the rest of the package builds on; this module imports no other.

``elements``, ``molecules`` and ``model`` all need the error type and the
diagnostic record, so they live here, below all three, and ``model`` can
import ``elements`` and ``molecules`` at the top like any other module.
"""

from importlib.metadata import version

#: Partial specific volume of a typical protein, cm³/g.
PROTEIN_SPECIFIC_VOLUME = 0.73

#: The severities a ``Diagnostic`` can have, mildest first.
SEVERITIES = ('info', 'warning')


class StructureError(ValueError):
    """A structure file or array set cannot be used for a calculation."""


def package_version(name: str) -> str | None:
    """The installed version of the distribution *name*, or None."""
    try:
        return version(name)
    except Exception:  # noqa: BLE001 - a missing distribution must not block a settings record
        return None


def residue_names(names, option: str) -> list[str]:
    """One residue name or an iterable of them, as a list.

    Raises:
        ValueError: If *names* is neither, naming the *option* it was given for.
    """
    if isinstance(names, str):
        return [names]
    try:
        values = list(names)
    except TypeError:
        values = None
    if values is None or not all(isinstance(name, str) for name in values):
        raise ValueError(f'{option} must be residue names, got {names!r}.')
    return [str(name) for name in values]


class Diagnostic(str):
    """One thing a reader or a calculation reports: its text, severity and origin.

    The ``notes`` of a ``Structure`` and a ``ScatteringProfile`` are these
    records, and every warning or log line the package prints about them is
    made from the record by one function. A diagnostic *is* its text as a
    ``str``, so notes still read, join, compare and serialise as plain strings;
    equality and hashing use the text alone.

    Attributes:
        severity: ``'warning'`` when the result may not be what was meant and
            the user can act on it; it is raised as a ``UserWarning``.
            ``'info'`` for a documented policy, applied as described; it is
            logged.
        origin: The step that raised it: ``'read_pdb'``, ``'read_gro'``,
            ``'read_sld'``, ``'read_vtk'``, ``'load'``, ``'calculate'``,
            ``'set_params'`` or ``'from_settings'``; ``'given'`` for a plain
            string passed in as a note.
    """

    severity: str
    origin: str

    def __new__(cls, text: str, severity: str = 'warning', origin: str = 'given') -> 'Diagnostic':
        if severity not in SEVERITIES:
            raise ValueError(f'severity must be one of {", ".join(SEVERITIES)}, got {severity!r}.')
        self = super().__new__(cls, text)
        self.__dict__.update(severity=severity, origin=str(origin))
        return self

    @classmethod
    def of(cls, note) -> 'Diagnostic':
        """*note* as a diagnostic: itself, or a plain string as a ``'given'`` warning."""
        return note if isinstance(note, cls) else cls(str(note))

    @property
    def text(self) -> str:
        """The text as a plain ``str``."""
        return str.__str__(self)

    def __setattr__(self, name: str, value) -> None:
        raise AttributeError('A Diagnostic cannot be changed.')

    def __delattr__(self, name: str) -> None:
        raise AttributeError('A Diagnostic cannot be changed.')

    def __reduce__(self):
        return type(self), (self.text, self.severity, self.origin)

    def __repr__(self) -> str:
        return f'Diagnostic({self.text!r}, severity={self.severity!r}, origin={self.origin!r})'
