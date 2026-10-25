"""Exception classes for expected EviTEA failures."""


class EviTEAError(Exception):
    """Base class for expected EviTEA failures."""


class InputContractError(EviTEAError, ValueError):
    """A user-provided input violates its declared format or semantics."""


class OutputStateError(InputContractError):
    """An output path is unsafe or incompatible with the requested operation."""


class PairSynchronizationError(InputContractError):
    """Paired reads cannot be consumed in strict lockstep."""


class RawInputValidationError(InputContractError):
    """Raw biological input is malformed or unsupported."""


class AssetCompatibilityError(EviTEAError):
    """Prepared assets do not share the required identities."""


class BackendCompatibilityError(EviTEAError):
    """A backend does not match its approved semantic profile."""


class ArchitectureInvariantError(EviTEAError, ValueError):
    """Internal state breaks a rule of the program design."""


class ResumeStateError(EviTEAError):
    """Checkpoint state is incomplete or incompatible."""


class RegressionCompatibilityError(EviTEAError):
    """A regression adapter cannot preserve the required identity."""


class OperationalResourceError(EviTEAError):
    """A process, storage, memory or external-resource operation failed."""


class OutputFinalizationError(EviTEAError):
    """Writing the final outputs and the provenance record failed."""
