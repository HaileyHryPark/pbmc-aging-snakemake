rule convert_cima_atac_h5ad_to_seu:
	input:
		data="data/cima_atac_data/cima_atac.h5ad",
	output:
		seu="data/cima_atac_data/cima_atac.rds",
	conda: "../env/anndata.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/cima_atac_data/convert_h5ad_to_seu.R"

rule extract_and_plot_metadata_cima_atac:
	input:
		seu="data/cima_atac_data/cima_atac.rds",
	output:
		meta="data/cima_atac_data/cima_atac_metadata_prefilter.csv",
		plot1="plots/cima_atac_data/cima_atac_qcplot1.svg",
		plot2="plots/cima_atac_data/cima_atac_qcplot2.svg",
		plot3="plots/cima_atac_data/cima_atac_qcplot3.svg",
	#conda: "../env/cima_atac_data.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_atac_data/extract_and_plot_metadata.R"

rule filter_and_normalize_data_cima_atac:
	input:
		seu="data/cima_atac_data/cima_atac.rds",
		meta="data/cima_atac_data/cima_atac_metadata_prefilter.csv",
	output:
		res="data/cima_atac_data/cima_atac_filtered.h5Seurat",
		meta="data/cima_atac_data/cima_atac_metadata_postfilter.csv",
	#conda: "../env/cima_atac_data.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_atac_data/filter_and_normalize_data.R"

