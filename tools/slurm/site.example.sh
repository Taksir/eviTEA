# Site settings (copy to $EVITEA_OUT/site.sh and edit). Paths are absolute or built from EVITEA_WORK.
SIF=$EVITEA_WORK/container/evitea.sif                 # Apptainer image with python, mappy, parasail, numpy
MSBWT_PATH=                                           # folder with the msbwt python package, if it is not in the image
REFFA=$EVITEA_WORK/reference/genome.fa                # reference genome (its checksum goes in the job log)
ASSETS=$EVITEA_WORK/reference/assets                  # prepared reference assets (evitea prepare-reference output)
PREVIOUS_SRC=$EVITEA_WORK/evitea_previous/src         # src folder of the previous version (a checkout of the repository's checkpoint branch)
SAMPLES=$EVITEA_OUT/samples.tsv                       # tab-separated: name, MSBWT index, library run folder, contigs
# samples.tsv columns 2 and 3 are relative to EVITEA_WORK, for example:
#   sampleA<TAB>data/sampleA/msbwt<TAB>runs/sampleA/call_out<TAB>chr2L,chr2R,chr3L,chr3R,chrX
