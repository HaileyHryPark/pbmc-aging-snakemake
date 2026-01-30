rule extract_metadata_cima_rna:
	input:
		"data/cima_rna_data_prep/cima_rna.h5ad",
	output:
		"data/cima_rna_data_prep/cima_rna_qced.h5ad",
		"data/cima_rna_data_prep/cima_rna_metadata.csv",
		"data/cima_rna_data_prep/cima_rna_ensembl_to_symbol.csv",
	params: dataset="cima_rna"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/extract_metadata.py"

rule write_initial_metadata_table_cima_rna:
	input:
		meta="data/cima_rna_data_prep/cima_rna_metadata.csv",
	output:
		table="tables/cima_rna_data_prep/cima_rna_data_prep_initial_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "04:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/write_metadata_table.R"

rule filter_data_cima_rna:
	input:
		"data/cima_rna_data_prep/cima_rna_qced.h5ad",
		"tables/cima_rna_data_prep/cima_rna_data_prep_initial_metadata_table.csv",
	output:
		"data/cima_rna_data_prep/cima_rna_filtered.h5ad",
		"data/cima_rna_data_prep/cima_rna_filtered_others.h5ad",
		"data/cima_rna_data_prep/cima_rna_metadata_filtered.csv",
	params: dataset="cima_rna"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "40:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/filter_data.py"

rule write_final_metadata_table_cima_rna:
	input:
		meta="data/cima_rna_data_prep/cima_rna_metadata_filtered.csv",
	output:
		table="tables/cima_rna_data_prep/internal_data_final_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/write_metadata_table.R"

rule split_h5ad_by_donor_cima_rna:
	input:
		"data/cima_rna_data_prep/cima_rna_filtered.h5ad",
	output:
		"data/cima_rna_data_prep/split_h5ad_by_donor_log_cima_rna.txt"
	params: dataset="cima_rna"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/split_h5ad_by_donor.py"

# manual_seu_conversion.R
rule convert_h5ad_to_seu_cima_rna:
	input:
		data="data/cima_rna_data_prep/cima_rna_filtered_{split}.h5ad",
	output:
		seu="data/cima_rna_data_prep/cima_rna_filtered_{split}.rds",
	conda: "../env/anndata.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/convert_h5ad_to_seu.R"

## Had to install azimuth manually and use server singularity (server outdated - old glibc); --use-singularity
rule run_azimuth_by_split_cima_rna:
	input:
		data="data/cima_rna_data_prep/cima_rna_filtered_{split}.rds",
	output:
		annotated="data/cima_rna_data_prep/cima_rna_{split}_azimuth_annotated.rds",
		plot="plots/cima_rna_data_prep/cima_rna_{split}_azimuth_annotation.pdf",
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/run_azimuth_by_split.R"

rule subset_donors_by_main_celltypes_cima_rna:
	input:
		data="data/cima_rna_data_prep/cima_rna_{split}_azimuth_annotated.rds",
	output:
		subset="data/cima_rna_data_prep/cima_rna_{split}_processed.rds",
		log="tables/cima_rna_data_prep/cima_rna_{split}_excluded_donors_log.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/subset_donors_by_main_celltypes.R"

rule summarize_final_data_included_cima_rna:
	input:
		data=expand("data/cima_rna_data_prep/cima_rna_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 10)]),
	output:
		summary="tables/cima_rna_data_prep/final_data_included_summary.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/cima_rna_data_prep/summarize_final_data_included.R"
		

