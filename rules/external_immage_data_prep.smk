rule extract_metadata_ext_immage:
	input:
		"data/external_immage_data_prep/immage.h5ad",
	output:
		"data/external_immage_data_prep/immage_qced.h5ad",
		"data/external_immage_data_prep/immage_metadata.csv",
		"data/external_immage_data_prep/immage_ensembl_to_symbol.csv",
	params: dataset="immage"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/extract_metadata.py"

rule write_initial_metadata_table_ext_immage:
	input:
		data="data/external_immage_data_prep/immage_metadata.csv",
	output:
		table="tables/external_immage_data_prep/external_data_initial_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/write_metadata_table.R"

rule filter_data_ext_immage:
	input:
		"data/external_immage_data_prep/immage_qced.h5ad",
		"tables/external_immage_data_prep/external_data_initial_metadata_table.csv",
	output:
		"data/external_immage_data_prep/immage_filtered.h5ad",
		"data/external_immage_data_prep/immage_metadata_filtered.csv",
	params: dataset="immage"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 400, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/filter_data.py"

rule write_final_metadata_table_ext_immage:
	input:
		data="data/external_immage_data_prep/immage_metadata_filtered.csv",
	output:
		table="tables/external_immage_data_prep/external_data_final_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/write_metadata_table.R"

rule split_h5ad_by_sample_ext_immage:
	input:
		"data/external_immage_data_prep/immage_filtered.h5ad",
	output:
		"data/external_immage_data_prep/split_h5ad_by_sample_log_immage.txt"
	params: dataset="immage"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/split_h5ad_by_sample.py"

# manual_seu_conversion.R
rule convert_h5ad_to_seu_ext_immage:
	input:
		data="data/external_immage_data_prep/immage_filtered_{split}.h5ad",
	output:
		seu="data/external_immage_data_prep/immage_filtered_{split}.rds",
	conda: "../env/anndata.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/convert_h5ad_to_seu.R"

## Had to install manually and use server singularity (server outdated - old glibc); --use-singularity
rule run_azimuth_by_split_ext_immage:
	input:
		data="data/external_immage_data_prep/immage_filtered_{split}.rds",
	output:
		annotated="data/external_immage_data_prep/immage_{split}_azimuth_annotated.rds",
		plot="plots/external_immage_data_prep/immage_{split}_azimuth_annotation.pdf",
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/run_azimuth_by_split.R"

rule subset_samples_by_main_celltypes_ext_immage:
	input:
		data="data/external_immage_data_prep/immage_{split}_azimuth_annotated.rds",
	output:
		subset="data/external_immage_data_prep/immage_{split}_processed.rds",
		log="tables/external_immage_data_prep/immage_{split}_excluded_samples_log.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/subset_samples_by_main_celltypes.R"

rule summarize_final_data_included_ext_immage:
	input:
		data=expand("data/external_immage_data_prep/immage_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 25)]),
	output:
		summary="tables/external_immage_data_prep/final_data_included_summary.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/summarize_final_data_included.R"
		
rule get_full5ct_pseudobulk_mat_immage:
	input:
		data="data/external_immage_data_prep/immage_{split}_processed.rds",
	output:
		pb="data/external_immage_data_prep/immage_{split}_full5ct_pseudobulk_data.csv",
	params: dataset="immage"
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/get_full5ct_pseudobulk_mat.R"

rule merge_full5ct_pseudobulk_data_immage:
	input:
		expand("data/external_immage_data_prep/immage_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 25)]),
	output:
		"data/external_immage_data_prep/full5ct_pseudobulk_data_immage.csv",
		"tables/external_immage_data_prep/full5ct_pseudobulk_data_immage_column_summary.txt",
		"plots/external_immage_data_prep/full5ct_pseudobulk_data_immage_pca.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/merge_full5ct_pseudobulk_data.py"

rule subset_allexp5ct_pseudobulk_immage:
	input:
		"data/external_immage_data_prep/full5ct_pseudobulk_data_immage.csv",
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv"
	output:
		"data/external_immage_data_prep/allexp5ct_pseudobulk_data_immage.csv",
		"data/external_immage_data_prep/allexp5ct_pseudobulk_data_immage_combined.csv",
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_data_prep/subset_allexp5ct_pseudobulk.py"

rule plot_sample_distribution_immage:
        input:
                data="data/external_immage_data_prep/full5ct_pseudobulk_data_immage.csv",
        output:
                plot="plots/external_immage_data_prep/immage_sample_distribution.svg",
        conda: "../env/internal_downstream.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 40, walltime = "05:00:00", queue = "super"
        script:
                "../scripts/external_immage_data_prep/plot_sample_distribution.R"

