# Documentation

The clean package currently exposes four command groups:

1. `prepare-reference` validates a host FASTA, optional repeat annotation, and
   binds a Bowtie2 index to a prepared-reference manifest.
2. `prepare-te` validates the TE library/family map, rebuilds the safe terminal
   and INTERNAL21 languages from the hash-bound frozen G1A candidate universe,
   and binds them to one prepared reference.
3. `prepare-sample` validates synchronized paired FASTQ, assigns ordinal-zero
   physical-pair identities, and builds a checked msBWT sample asset.
4. `validate` verifies install, prepared assets, or finalized run artifacts
   without truth data or scoring.

The default end-to-end engine is available as `evitea call`: it runs the clean
incumbent caller followed by the four evidence-extension layers on one prepared
sample. See [`CLI_AND_DATA.md`](CLI_AND_DATA.md) for inputs, commands and
outputs, and [`../validation/EVITEA_V1_VALIDATION_STATUS.md`](../validation/EVITEA_V1_VALIDATION_STATUS.md)
for what has been validated.

Historical research-stage names appear only in planning and validation
provenance. The production naming map is in [`NAMING.md`](NAMING.md).
