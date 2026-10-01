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

Run the clean incumbent (R9E) engine end to end on one development sample and
compare every stage with the frozen Route-B artifacts (resumable; needs the author's private development
datasets, see `validation/README.md`):

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_INCUMBENT.ps1 -Sample ISO1 -Work <work-dir>
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_INCUMBENT.ps1 -Sample AKA017 -Work <work-dir>
```

`-Sample ISO1` (default) compares with `TE_locator_route_b_v2\reports` and the
frozen ISO1 R9E (192 rows, `b3e48681…`) and uses the ISO1 recall-repair query
names `{pid}`; work directory `<work-dir>`.
Pass `-Workers N` and any `run_incumbent_validation.py` option after the named
parameters.  Verdict: `INCUMBENT_EQUIVALENCE_PASS` in
`<Work>\INCUMBENT_VALIDATION_REPORT.json`.  If native Windows Bowtie2 rejects
the all-numeric ISO1 recall-repair names, pass `--bowtie2 <WSL/Linux bowtie2 wrapper>`.

Check the v3 language builders (codebooks) against the frozen tables (~10 min,
no sample data):

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_LANGUAGES.ps1 -Cache <cache-dir> -RouteBV3 <route-b-v3-dir>
```

Verdict: `LANGUAGE_BUILDER_EQUIVALENCE_PASS` in
`<cache-dir>\LANGUAGE_EQUIVALENCE.json`.

## End-to-end `evitea call` validation (both samples)

Runs the public pipeline exactly as a user would: `prepare-reference`
(adopting the frozen dm6 Bowtie2 index), `prepare-te --build-languages`
(building both v3 codebooks with EviTEA and checking them against the frozen
codebooks; `-Languages supplied` binds the frozen files instead), then `evitea call --msbwt-index <msBWT> --profile <sample>` on
ISO1 and AKA017, and checks the SHA-256 after every layer (R9E incumbent,
bilateral, asymmetric, full memory, INTERNAL21) plus `calls.native.tsv`
against the frozen callsets (AKA017 838 `6432c8e1…`, ISO1 194 `9e538d74…`).
Resumable; rerun the same command after an interruption.

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_EVITEA_CALL.ps1 -Work <work-dir> -RouteBV3 <route-b-v3-dir>
```

Options: `-Samples ISO1` (one sample), `-Work <work-dir>`,
`-Threads 8`, `-Workers N`, `-EdlibPython <python with edlib>`.  Verdict file:
`<Work>\EVITEA_CALL_END_TO_END_EQUIVALENCE.json`
(`EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS`).

## Paths you must set yourself

These validation launchers contain no machine-specific paths. `-Work`, `-Cache`,
`-RouteBV2`/`-RouteBV3` (and `-EdlibPython` for the probe) are arguments you pass, and the
underlying Python scripts read `EVITEA_DATA` (the folder holding the development
datasets) from the environment. See `validation/README.md`.
