"""Adapter for the ``pymsbwt`` public query API (module ``MUSCython.MultiStringBWTCython``)."""

from __future__ import annotations

import hashlib
import logging
import os
import re
import shutil
import stat
from collections.abc import Iterable
from importlib import metadata
from pathlib import Path
from typing import Any, Callable

from ..contracts import ReadMate
from ..errors import (
    AssetCompatibilityError,
    BackendCompatibilityError,
    InputContractError,
)
from ..identity import iter_physical_pairs
from ..inputs import iter_paired_fastq
from .exact_query import (
    PreparedSample,
    ReadIdentityMap,
    make_read_identity_map,
    prepare_sample_manifest,
    write_read_identity_map,
)
from .protocols import MapperRead, RawReadHit, ReadIdentity

MSBWT_DISTRIBUTION = "pymsbwt"
# Lowest pymsbwt release accepted; every release builds the same msBWT byte for byte.
MSBWT_MIN_VERSION = (0, 5)
# The earlier distribution of the same module, accepted as a fallback at its validated version.
LEGACY_MSBWT_DISTRIBUTION = "msbwt-modern3"
LEGACY_MSBWT_VERSION = "0.3.0"
INDEX_DIRECTORY = "index"
IDENTITY_MAP_NAME = "read_identity_map.json"
REQUIRED_INDEX_FILES = ("about.npy", "fmIndex.npy", "msbwt.npy", "totalCounts.npy")
_MSBWT_ALPHABET = frozenset("ACGTN")


def _make_tree_writable(path: Path) -> None:
    """Clear Windows read-only flags left by the compiled builder."""

    if os.name != "nt" or not path.exists():
        return
    for target in (path, *path.rglob("*")):
        try:
            os.chmod(target, target.stat().st_mode | stat.S_IWRITE)
        except FileNotFoundError:
            continue


def _remove_builder_staging(path: Path) -> None:
    """Remove Windows builder scratch files, including read-only artifacts."""

    if not path.exists():
        return

    def make_writable_and_retry(
        function: Callable[..., Any], target: str, _: BaseException
    ) -> None:
        os.chmod(target, stat.S_IWRITE)
        function(target)

    shutil.rmtree(path, onexc=make_writable_and_retry)


def _normalize_query(sequence: str) -> str:
    if not isinstance(sequence, str) or not sequence:
        raise InputContractError("exact query must be non-empty text")
    normalized = sequence.upper()
    if any(base not in "ACGT" for base in normalized):
        raise InputContractError("exact query must contain only A/C/G/T bases")
    return normalized


def load_msbwt(index_dir: str | Path, *, use_memmap: bool = True) -> Any:
    """Load an index lazily so importing EviTEA does not require msBWT."""

    try:
        from MUSCython.MultiStringBWTCython import loadBWT
    except ImportError as exc:
        raise BackendCompatibilityError(
            "pymsbwt is unavailable; install it with `pip install pymsbwt`"
        ) from exc
    try:
        return loadBWT(str(Path(index_dir)), useMemmap=use_memmap)
    except Exception as exc:
        raise BackendCompatibilityError(
            f"failed to load msBWT index: {index_dir}"
        ) from exc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _release_tuple(version: str) -> tuple[int, ...]:
    """Leading numeric release parts: ``"0.5.0"`` -> (0, 5, 0), ``"0.6rc1"`` -> (0, 6)."""

    match = re.match(r"\d+(?:\.\d+)*", version)
    return tuple(int(part) for part in match.group().split(".")) if match else ()


def installed_msbwt_distribution() -> tuple[str, str]:
    """(distribution, version) of the installed msBWT package: pymsbwt, else the legacy package."""

    for name in (MSBWT_DISTRIBUTION, LEGACY_MSBWT_DISTRIBUTION):
        try:
            return name, metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    raise BackendCompatibilityError("pymsbwt is unavailable; install it with `pip install pymsbwt`")


def _installed_uniform_builder() -> tuple[Callable[..., None], str, str, Path]:
    try:
        distribution, installed_version = installed_msbwt_distribution()
        from MUSCython import MultiStringBWTCython
    except ImportError as exc:
        raise BackendCompatibilityError(
            "pymsbwt with compiled builders is unavailable"
        ) from exc
    if distribution == MSBWT_DISTRIBUTION:
        if _release_tuple(installed_version) < MSBWT_MIN_VERSION:
            raise BackendCompatibilityError(
                f"pymsbwt {installed_version} is older than the required "
                f"{'.'.join(map(str, MSBWT_MIN_VERSION))}"
            )
    elif installed_version != LEGACY_MSBWT_VERSION:
        raise BackendCompatibilityError(
            f"{LEGACY_MSBWT_DISTRIBUTION} {installed_version} is not the validated "
            f"{LEGACY_MSBWT_VERSION}; install pymsbwt instead"
        )
    builder = getattr(MultiStringBWTCython, "createMSBWTFromFastq", None)
    module_path = Path(getattr(MultiStringBWTCython, "__file__", ""))
    if not callable(builder) or not module_path.is_file():
        raise BackendCompatibilityError(
            f"{distribution} lacks the compiled uniform FASTQ builder"
        )
    return builder, distribution, installed_version, module_path


def _identity_map_from_about(
    about_path: Path,
    *,
    pair_count: int,
) -> ReadIdentityMap:
    try:
        import numpy as np

        about = np.load(about_path, mmap_mode="r")
    except (ImportError, OSError, ValueError) as exc:
        raise AssetCompatibilityError(
            "msBWT about.npy is missing or unreadable"
        ) from exc
    names = about.dtype.names
    if about.ndim == 1 and names and len(names) >= 2:
        fields = names[:2]
        rows = ((int(row[fields[0]]), int(row[fields[1]])) for row in about)
    elif about.ndim == 2 and about.shape[1] >= 2:
        rows = ((int(row[0]), int(row[1])) for row in about)
    else:
        raise AssetCompatibilityError("msBWT about.npy has an unsupported schema")
    identities = []
    for backend_read_id, (file_id, read_ordinal0) in enumerate(rows):
        if file_id not in (0, 1) or not 0 <= read_ordinal0 < pair_count:
            raise AssetCompatibilityError(
                "msBWT about.npy contains invalid read identity"
            )
        identities.append(
            ReadIdentity(
                backend_read_id=backend_read_id,
                pair_id=read_ordinal0,
                mate=ReadMate.READ1 if file_id == 0 else ReadMate.READ2,
            )
        )
    identity_map = make_read_identity_map(identities)
    observed = {(row.pair_id, row.mate) for row in identity_map.identities.values()}
    expected = {
        (pair_id, mate)
        for pair_id in range(pair_count)
        for mate in (ReadMate.READ1, ReadMate.READ2)
    }
    if observed != expected or len(identity_map.identities) != pair_count * 2:
        raise AssetCompatibilityError(
            "msBWT about.npy does not preserve every physical pair and mate"
        )
    return identity_map


def build_prepared_msbwt_sample(
    read1: str | Path,
    read2: str | Path,
    output_dir: str | Path,
    *,
    sample_id: str,
    threads: int,
    builder: Callable[..., None] | None = None,
    initializer: Callable[[Path], Any] | None = None,
    backend_version: str | None = None,
    builder_module_path: str | Path | None = None,
) -> PreparedSample:
    """Build the validated uniform-FASTQ profile and bind dollar IDs to pairs."""

    if not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        raise InputContractError("msBWT builder threads must be a positive integer")
    root = Path(output_dir).expanduser().resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise InputContractError("prepared-sample output must be absent or empty")
    root.mkdir(parents=True, exist_ok=True)

    if builder is None:
        builder, distribution, installed_version, installed_module_path = _installed_uniform_builder()
    else:
        distribution = MSBWT_DISTRIBUTION
        installed_version = backend_version or "test-double"
        installed_module_path = (
            Path(builder_module_path) if builder_module_path else None
        )

    staging = root / ".builder-input"
    staging.mkdir()
    normalized_read1 = staging / "read1.fastq"
    normalized_read2 = staging / "read2.fastq"
    pair_count = 0
    uniform_length: int | None = None
    try:
        with (
            normalized_read1.open("w", encoding="ascii", newline="\n") as out1,
            normalized_read2.open("w", encoding="ascii", newline="\n") as out2,
        ):
            pairs = iter_physical_pairs(
                iter_paired_fastq(read1, read2), sample_id=sample_id
            )
            for pair in pairs:
                for read, output, mate in (
                    (pair.read1, out1, 1),
                    (pair.read2, out2, 2),
                ):
                    if set(read.sequence) - _MSBWT_ALPHABET:
                        raise InputContractError(
                            "msBWT input supports only A/C/G/T/N FASTQ bases"
                        )
                    if uniform_length is None:
                        uniform_length = len(read.sequence)
                    elif len(read.sequence) != uniform_length:
                        raise BackendCompatibilityError(
                            "validated msBWT identity profile requires "
                            "uniform read length"
                        )
                    output.write(
                        f"@pair_{pair.pair_id}/{mate}\n{read.sequence}\n+\n"
                        f"{'I' * len(read.sequence)}\n"
                    )
                pair_count += 1
        if pair_count == 0:
            raise InputContractError("paired FASTQ input is empty")

        index_dir = root / INDEX_DIRECTORY
        index_dir.mkdir()
        try:
            builder(
                [str(normalized_read1), str(normalized_read2)],
                str(index_dir),
                threads,
                True,
                logging.getLogger("evitea.msbwt"),
            )
        except Exception as exc:
            raise BackendCompatibilityError("msBWT uniform FASTQ build failed") from exc

        _make_tree_writable(root)
        try:
            initialized_index = (initializer or load_msbwt)(index_dir)
            if initialized_index is None:
                raise BackendCompatibilityError(
                    "msBWT index initialization returned no index"
                )
            del initialized_index
        except BackendCompatibilityError:
            raise
        except Exception as exc:
            raise BackendCompatibilityError(
                "msBWT index initialization failed"
            ) from exc
    finally:
        _make_tree_writable(root)
        _remove_builder_staging(staging)

    identity_map = _identity_map_from_about(
        root / INDEX_DIRECTORY / "about.npy", pair_count=pair_count
    )
    write_read_identity_map(root / IDENTITY_MAP_NAME, identity_map)
    index_files = tuple(
        (Path(INDEX_DIRECTORY) / filename).as_posix()
        for filename in REQUIRED_INDEX_FILES
    )
    missing_index_files = [
        relative_path
        for relative_path in index_files
        if not (root / relative_path).is_file()
    ]
    if missing_index_files:
        raise AssetCompatibilityError(
            "msBWT builder omitted required index files: "
            + ", ".join(missing_index_files)
        )
    module_sha256 = (
        _sha256(installed_module_path)
        if installed_module_path is not None and installed_module_path.is_file()
        else None
    )
    return prepare_sample_manifest(
        root,
        sample_id=sample_id,
        pair_count=pair_count,
        read1_path=read1,
        read2_path=read2,
        identity_map_relative_path=IDENTITY_MAP_NAME,
        backend_profile={
            "distribution": distribution,
            "version": installed_version,
            "builder": "MUSCython.MultiStringBWTCython.createMSBWTFromFastq",
            "builder_module_sha256": module_sha256,
            "uniform_length": uniform_length,
            "input_file_order": ["read1", "read2"],
            "identity_source": f"{INDEX_DIRECTORY}/about.npy",
            "index_directory": INDEX_DIRECTORY,
            "threads": threads,
        },
        index_relative_paths=index_files,
    )


class MsbwtExactQueryBackend:
    """Exact-query observations backed by a prebuilt msBWT index.

    ``getSequenceDollarID`` is used only as the backend read key. Pair and mate
    identity come exclusively from the prepared explicit identity map.
    """

    def __init__(
        self,
        index_dir: str | Path,
        identity_map: ReadIdentityMap,
        *,
        native_bwt: Any | None = None,
        loader: Callable[[str | Path], Any] | None = None,
    ) -> None:
        self.index_dir = Path(index_dir)
        self.identity_map = identity_map
        if native_bwt is not None and loader is not None:
            raise InputContractError("provide native_bwt or loader, not both")
        self._bwt = (
            native_bwt
            if native_bwt is not None
            else (loader or load_msbwt)(self.index_dir)
        )
        for method in ("findIndicesOfStr", "getSequenceDollarID"):
            if not callable(getattr(self._bwt, method, None)):
                raise BackendCompatibilityError(
                    f"msBWT object lacks required public method {method}"
                )

    def build(self, reads: Iterable[MapperRead], output_dir: str | Path) -> None:
        del reads, output_dir
        raise BackendCompatibilityError(
            "msBWT building requires a separately validated builder profile"
        )

    def _interval(self, sequence: str) -> tuple[str, int, int]:
        query = _normalize_query(sequence)
        try:
            lo, hi = self._bwt.findIndicesOfStr(query)
            lo_int, hi_int = int(lo), int(hi)
        except Exception as exc:
            raise BackendCompatibilityError("msBWT exact query failed") from exc
        if lo_int < 0 or hi_int < lo_int:
            raise BackendCompatibilityError(
                "msBWT returned an invalid half-open interval"
            )
        return query, lo_int, hi_int

    def count_exact(self, sequence: str) -> int:
        _, lo, hi = self._interval(sequence)
        return hi - lo

    def locate_reads(self, sequence: str) -> Iterable[RawReadHit]:
        query, lo, hi = self._interval(sequence)
        for occurrence_ordinal0, row in enumerate(range(lo, hi)):
            try:
                backend_read_id = int(self._bwt.getSequenceDollarID(row))
            except Exception as exc:
                raise BackendCompatibilityError(
                    f"msBWT could not recover read identity for row {row}"
                ) from exc
            identity = self.identity_map.resolve(backend_read_id)
            yield RawReadHit(
                sequence=query,
                backend_row=row,
                backend_read_id=backend_read_id,
                pair_id=identity.pair_id,
                mate=identity.mate,
                occurrence_ordinal0=occurrence_ordinal0,
            )

    def recover_read(self, backend_read_id: int) -> str:
        """Recover one immutable input read from its bound msBWT dollar ID."""

        self.identity_map.resolve(backend_read_id)
        recover = getattr(self._bwt, "recoverString", None)
        if not callable(recover):
            raise BackendCompatibilityError(
                "msBWT object lacks required public method recoverString"
            )
        try:
            recovered = recover(backend_read_id)
            if isinstance(recovered, bytes):
                value = recovered.decode("ascii")
            elif isinstance(recovered, str):
                value = recovered
            else:
                raise TypeError("recoverString returned non-text data")
        except Exception as exc:
            raise BackendCompatibilityError(
                f"msBWT could not recover read {backend_read_id}"
            ) from exc
        if not value.startswith("$") or value.count("$") != 1:
            raise BackendCompatibilityError(
                "msBWT recovered read lacks exactly one leading sentinel"
            )
        sequence = value[1:].upper()
        if not sequence or any(base not in "ACGTN" for base in sequence):
            raise BackendCompatibilityError(
                "msBWT recovered read is empty or has invalid bases"
            )
        return sequence

    def recover_sample_read_index(self, sample_id: str):
        """Build the physical-read projection index without raw FASTQ sidecars."""

        from ..evidence.projection import SampleReadIndex

        sequences = {}
        for backend_read_id, identity in self.identity_map.identities.items():
            key = (identity.pair_id, identity.mate)
            if key in sequences:
                raise BackendCompatibilityError(
                    "msBWT identity map contains duplicate physical-read identity"
                )
            sequences[key] = self.recover_read(backend_read_id)
        return SampleReadIndex(sample_id, sequences)
