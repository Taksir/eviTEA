# Windows helpers

`RUN_EVITEA.ps1` invokes the same package CLI and contains no scientific logic.
For example:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\RUN_EVITEA.ps1 version --json
```

Pass `-Python C:\path\to\python.exe` before the remaining EviTEA arguments when
the supported interpreter is not the default `python` on `PATH`.

Run the source regression suite and wheel/sdist build on an ordinary local
Windows shell with:

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_LOCAL.ps1
```

Add `-RunLiveMsbwt` to include the compiled msBWT integration test. The script
does not alter package or scientific configuration; it exists because the
managed development environment may deny Python build hooks access to their
private temporary directories.

Run the library layer end to end on one development sample and
compare every stage with the fixed reference artifacts (resumable; requires the author's private development
datasets, see `validation/README.md`):

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_INCUMBENT.ps1 -Sample ISO1 -Work <work-dir>
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_INCUMBENT.ps1 -Sample AKA017 -Work <work-dir>
```

`-Sample ISO1` (default) compares with `TE_locator_route_b_v2\reports` and the
fixed ISO1 library-layer result (192 rows, `b3e48681...`) and uses the ISO1 wide-anchor query
names `{pid}`; work directory `<work-dir>`.
Pass `-Workers N` and any `run_incumbent_validation.py` option after the named
parameters.  Result: `INCUMBENT_EQUIVALENCE_PASS` in
`<Work>\INCUMBENT_VALIDATION_REPORT.json`.  If native Windows Bowtie2 rejects
the all-numeric ISO1 wide-anchor names, pass `--bowtie2 <WSL/Linux bowtie2 wrapper>`.

Check the language builders (codebooks) against the fixed tables (~10 min,
no sample data):

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_LANGUAGES.ps1 -Cache <cache-dir> -RouteBV3 <route-b-v3-dir>
```

Result: `LANGUAGE_BUILDER_EQUIVALENCE_PASS` in
`<cache-dir>\LANGUAGE_EQUIVALENCE.json`.

## End-to-end `evitea call` validation (both samples)

Runs the public pipeline exactly as a user would: `prepare-reference`
(adopting the fixed dm6 Bowtie2 index), `prepare-te --build-languages`
(building both codebooks with EviTEA and checking them against the fixed
codebooks; `-Languages supplied` binds the fixed files instead), then `evitea call --msbwt-index <msBWT> --profile <sample>` on
ISO1 and AKA017, and checks the SHA-256 after every layer (library layer,
bilateral, asymmetric, full memory, internal 21-mer) plus `calls.native.tsv`
against the fixed call sets (AKA017 838 `6432c8e1...`, ISO1 194 `9e538d74...`).
Resumable; rerun the same command after an interruption.

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_EVITEA_CALL.ps1 -Work <work-dir> -RouteBV3 <route-b-v3-dir>
```

Options: `-Samples ISO1` (one sample), `-Work <work-dir>`,
`-Threads 8`, `-Workers N`, `-EdlibPython <python with edlib>`.  Result file:
`<Work>\EVITEA_CALL_END_TO_END_EQUIVALENCE.json`
(`EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS`).

## Paths you must set yourself

These validation launchers contain no machine-specific paths. `-Work`, `-Cache`,
`-RouteBV2`/`-RouteBV3` (and `-EdlibPython` for the probe) are arguments you pass, and the
underlying Python scripts read `EVITEA_DATA` (the folder holding the development
datasets) from the environment. See `validation/README.md`.
