# Validation assets

This tree is reserved for frozen regression assets, legacy adapters, and
optional quarantined scoring. Validation code may import `evitea`; production
code must not import this tree. No truth data, comparator outputs, or large
regression assets are included in the Phase 0 scaffold.

Reserved directories: `fixtures/`, `frozen_regression/`, `legacy_adapters/`,
`negative_regression/`, and `quarantined_scoring/`.

Implemented checks (validation-only; production never imports them):

* `incumbent_stage_checks/` — clean incumbent vs frozen stage tables
  (`run_incumbent_validation.py --sample AKA017|ISO1`; VM chunked ISO1 drivers
  `iso1_acquisition_vm.py`, `iso1_downstream_vm.py`).
* `v3_layer_checks/` — v3 extension layers and engines.
* `language_checks/` — v3 language (codebook) builders vs frozen tables.
* `end_to_end/` — the public `evitea call` pipeline, every layer hash.
* `frozen_regression/` — recorded layer hashes (oracle integrity).

## Paths you must set yourself

The validation scripts reproduce frozen results from the author's private
development datasets (two *Drosophila melanogaster* samples, AKA017 and ISO1)
and the frozen research outputs they were compared against. **Those datasets
and research outputs are not distributed with this repository**, and no
machine-specific path is hard-coded. If you have your own copy, put the paths
in yourself:

| What | How to set it |
|---|---|
| Folder with the sample data, TE library and reference (`EVITEA_DATA`) | environment variable, or `--data` |
| Frozen Route-B v2 / v3 research folders | `EVITEA_ROUTE_B_V2`, `EVITEA_ROUTE_B_V3`, or `--route-b-v2` / `--route-b-v3` |
| Work / cache directory | required `--work` / `--cache` argument (`-Work` / `-Cache` in the PowerShell launchers) |
| Python with `edlib` (only where noted) | `--edlib-python` / `-EdlibPython` |
| Path to this checkout, for the `incumbent_stage_checks` VM drivers | `EVITEA_REPO` |

PowerShell example:

```powershell
$env:EVITEA_DATA = "<path to your data folder>"
$env:EVITEA_ROUTE_B_V3 = "<path to the frozen route-b-v3 folder>"
.\tools\windows\VALIDATE_EVITEA_CALL.ps1 -Work <work-dir> -RouteBV3 $env:EVITEA_ROUTE_B_V3
```
