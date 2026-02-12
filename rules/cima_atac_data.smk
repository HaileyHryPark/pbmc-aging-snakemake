rule extract_and_plot_metadata_cima_atac:
	input:
		"data/cima_atac_data/cima_atac.h5ad",
	output:
		"data/cima_atac_data/cima_atac_metadata_prefilter.csv",
	conda: "../env/cima_atac_data.yaml"
	#singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_atac_data/extract_and_plot_metadata.py"

rule filter_and_normalize_data_cima_atac:
	input:
		"data/cima_atac_data/cima_atac.h5ad",
		"data/cima_atac_data/cima_atac_metadata_prefilter.csv",
	output:
		"data/cima_atac_data/cima_atac_processed.h5ad",
	conda: "../env/cima_atac_data.yaml"
	#singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_atac_data/filter_and_normalize_data.py"

