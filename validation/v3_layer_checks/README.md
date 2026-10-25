# Extension-layer checks

Validation only. Package code never imports anything here.

## `run_shadow_layer_equivalence.py`: integration equivalence

Feeds the recorded research inputs of the four integration steps (paired-end
anchor candidates, full-memory selected anchors, internal 21-mer locus record
and event instances) into `evitea.shadow_layers` on top of the library-layer
callset bytes. It requires byte identity at every layer plus exact agreement of
every integration decision and nearest-call attribution:

```text
AKA017  792 -> 801 -> 820 -> 835 -> 838   final 6432c8e1...f701
ISO1    192 -> 192 -> 193 -> 194 -> 194   final 9e538d74...92b0
```

```powershell
cd <path-to-EviTEA>
python validation\v3_layer_checks\run_shadow_layer_equivalence.py `
  --route-b-v2 "<route-b-v2-dir>" `
  --route-b-v3 "<route-b-v3-dir>"
```

Pass `--incumbent-aka017 <engine R9E_COMBINED_CALLS.tsv>` to build on the clean
library-layer engine output instead of the recorded library-layer file (the two
are byte-identical after `INCUMBENT_EQUIVALENCE_PASS`).

Scope: this isolates the integration semantics. It does **not** show that
EviTEA derives the layer inputs itself. The engine checks below cover that.

## `run_iso_procedure_on_aka017.py`: ISO1 procedure on AKA017 (laptop, hours)

Runs the fixed, self-contained ISO1 paired-end anchor transfer procedure
(imported by SHA-256) unchanged on AKA017. It compares its corrected clusters,
novel/nonresident targets, replicated exact-junction anchors, paired-end
support universe and candidates with the recorded AKA017 chain. Results:

* `ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_EXACTLY`
* `ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_UP_TO_CLUSTER_IDS`
* `ISO_PROCEDURE_DIVERGES_FROM_AKA017_FROZEN_CHAIN`

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\PROBE_ISO_PROCEDURE_AKA017.ps1 -Work <work-dir> -RouteBV3 <route-b-v3-dir> -EdlibPython <python-with-edlib>
```

Requires the edlib-capable interpreter used by the research scripts
(`-EdlibPython <path\to\python.exe>`, required) and an
msbwt-modern3 Python (`-Python`). The result is
`ISO_PROCEDURE_AKA017_COMPARISON.json` in the work directory.

## `run_v3_extension_equivalence.py`: clean extension engine (laptop)

Runs `evitea.extension_engine.run_extension` from the prepared msBWT of each
sample (AKA017 host read names `aka-historical`, ISO1 `iso`). It requires every
intermediate (paired-end anchor candidates; full-memory proposals, event record
and selected events; internal 21-mer witnesses, event instances and loci) to
equal the recorded tables, and the four layers to be byte-identical (AKA017
801/820/835/838, ISO1 192/193/194/194). The result is
`V3_EXTENSION_ENGINE_EQUIVALENCE_PASS` in `V3_EXTENSION_ENGINE_EQUIVALENCE.json`.

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_V3_EXTENSION.ps1 -Work <work-dir> -RouteBV2 <route-b-v2-dir> -RouteBV3 <route-b-v3-dir>
```
