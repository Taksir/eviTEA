# Documentation

- [`../README.md`](../README.md): what EviTEA does, installation, a quick start and what is validated.
- [`CLI_AND_DATA.md`](CLI_AND_DATA.md): every command, its inputs and outputs, and the values taken from the sample and the
  reference (read length, fragment size, primary contigs).
- [`REPEAT_ANNOTATION.md`](REPEAT_ANNOTATION.md): making the repeat annotation with RepeatMasker, and runs without a TE
  library (`build-repeat-map`, `evitea run --no-library`).
- [`LIBRARY_FREE_OPTIONS.md`](LIBRARY_FREE_OPTIONS.md): the switches of the library-free layer.
- [`NAMING.md`](NAMING.md): plain names for the research code's shorthand labels.
- [`../examples/README.md`](../examples/README.md): a worked example on synthetic data.

Commands, in the order of a run:

| command | what it does |
|---|---|
| `prepare-reference` | checks the reference FASTA and its repeat annotation, builds (or adopts) the Bowtie2 index |
| `prepare-te` | checks the TE library and family map, builds the TE codebooks (`--build-languages`) |
| `prepare-sample` | checks paired FASTQ and builds the sample's msBWT |
| `call` | the library calls |
| `library-free` | the library-free calls, added to the library calls (Linux/WSL) |
| `build-repeat-map`, `build-family-records`, `add-record-calls` | the steps of a run without a TE library |
| `run` | all of the above as one command, with a library or `--no-library` |
| `refine-breakpoints`, `call-tiers`, `filter-resident-copies` | optional views of a finished run |
| `validate`, `check-install`, `version` | checks prepared assets, finished runs and the installation |

What has been validated is in [`../validation/EVITEA_V1_VALIDATION_STATUS.md`](../validation/EVITEA_V1_VALIDATION_STATUS.md).
