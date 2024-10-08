nextflow run main.nf \
	--INDEX  /private/groups/kimlab/vikas/refs/metagenomics_index/ \
    --INPUT_FOLDER /private/groups/kimlab/vikas/metagenomics/publication/ncm_2024/input/time_tagged_input/CSF001/ \
	--LOW_COMPLEXITY_FILTER_NANOPORE \
	--ALIGN_ALL_FUNGI \
	--KRAKEN_PREFILTER \
    --PREFILTER_THRESHOLD 5 \
	--CHOPPER_MAXLENGTH 500000 \
	--CHOPPER_QUALITY 9 \
	--OUTPUT aPOMP_out_kraken \
	-with-docker ubuntu:18.04 \
	-with-report \
	-with-trace \
    -profile cluster \
    -resume