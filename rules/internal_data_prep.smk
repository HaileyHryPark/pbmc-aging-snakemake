rule prep_marina_data:
	input:
		"data/internal_data_prep/all_pbmcs/all_pbmcs_rna.h5ad",
		"data/internal_data_prep/all_pbmcs/all_pbmcs_metadata.csv",
	output:
		"data/internal_data_prep/marina.h5ad",
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/prep_marina_data.py"

rule extract_metadata:
	input:
		"data/internal_data_prep/{dataset}.h5ad",
	output:
		"data/internal_data_prep/{dataset}_qced.h5ad",
		"data/internal_data_prep/{dataset}_metadata.csv",
		"data/internal_data_prep/{dataset}_ensembl_to_symbol.csv",
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/extract_metadata.py"

rule write_initial_metadata_table:
	input:
		onek1k="data/internal_data_prep/onek1k_metadata.csv",
		aida="data/internal_data_prep/aida_metadata.csv",
		marina="data/internal_data_prep/marina_metadata.csv",
		perez="data/internal_data_prep/perez_metadata.csv",
	output:
		table="tables/internal_data_prep/internal_data_initial_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/write_metadata_table.R"

rule filter_data:
	input:
		"data/internal_data_prep/{dataset}_qced.h5ad",
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
	output:
		"data/internal_data_prep/{dataset}_filtered.h5ad",
		"data/internal_data_prep/{dataset}_filtered_others.h5ad",
		"data/internal_data_prep/{dataset}_metadata_filtered.csv",
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/filter_data.py"

rule write_final_metadata_table:
	input:
		onek1k="data/internal_data_prep/onek1k_metadata_filtered.csv",
		aida="data/internal_data_prep/aida_metadata_filtered.csv",
		marina="data/internal_data_prep/marina_metadata_filtered.csv",
		perez="data/internal_data_prep/perez_metadata_filtered.csv",
	output:
		table="tables/internal_data_prep/internal_data_final_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/write_metadata_table.R"

rule split_h5ad_by_donor:
	input:
		"data/internal_data_prep/{dataset}_filtered.h5ad",
	output:
		"data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt"
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/split_h5ad_by_donor.py"

# manual_seu_conversion.R
rule convert_h5ad_to_seu:
	input:
		data="data/internal_data_prep/{dataset}_filtered_{split}.h5ad",
	output:
		seu="data/internal_data_prep/{dataset}_filtered_{split}.rds",
	conda: "../env/anndata.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/convert_h5ad_to_seu.R"

## Had to install azimuth manually and use server singularity (server outdated - old glibc); --use-singularity
rule run_azimuth_by_split:
	input:
		data="data/internal_data_prep/{dataset}_filtered_{split}.rds",
	output:
		annotated="data/internal_data_prep/{dataset}_{split}_azimuth_annotated.rds",
		plot="plots/internal_data_prep/{dataset}_{split}_azimuth_annotation.pdf",
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/run_azimuth_by_split.R"

rule subset_donors_by_main_celltypes:
	input:
		data="data/internal_data_prep/{dataset}_{split}_azimuth_annotated.rds",
	output:
		subset="data/internal_data_prep/{dataset}_{split}_processed.rds",
		log="tables/internal_data_prep/{dataset}_{split}_excluded_donors_log.csv",
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/subset_donors_by_main_celltypes.R"

rule summarize_final_data_included:
	input:
		onek1k=expand("data/internal_data_prep/onek1k_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 11)]),
		aida=expand("data/internal_data_prep/aida_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 22)]),
		perez=expand("data/internal_data_prep/perez_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 5)]),
		marina=expand("data/internal_data_prep/marina_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 7)]),
	output:
		summary="tables/internal_data_prep/final_data_included_summary.csv",
		cellcount="tables/internal_data_prep/final_data_cellcount_df.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_data_prep/summarize_final_data_included.R"
		
rule plot_final_data_cellcount_per_dataset:
	input:
		cellcount="tables/internal_data_prep/final_data_cellcount_df.csv",
	output:
		plot="plots/internal_data_prep/final_data_cellcount_by_dataset.svg"	
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 20, walltime = "10:00:00", queue = "short"
	script:
		"../scripts/internal_data_prep/plot_final_data_cellcount_per_dataset.R"
		
			
