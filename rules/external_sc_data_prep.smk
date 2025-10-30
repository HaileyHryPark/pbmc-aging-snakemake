rule create_seurat_obj_sc:
	input:
		mtx="data/external_sc_data_prep/01.UMI.txt",
		metadata="data/external_sc_data_prep/03.Cell.Barcodes.txt",
	output:
		seu="data/external_sc_data_prep/sc_qced.rds",
		meta="data/external_sc_data_prep/sc_metadata.csv",
		annot="data/external_sc_data_prep/sc_annot.csv",
	conda: "../env/external_sc_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/create_seurat_obj_sc.R"

rule write_initial_metadata_table_ext_sc:
	input:
		meta="data/external_sc_data_prep/sc_metadata.csv",
	output:
		table="tables/external_sc_data_prep/external_data_initial_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/write_metadata_table.R"

rule filter_data_ext_sc:
	input:
		seu="data/external_sc_data_prep/sc_qced.rds",
		table="tables/external_sc_data_prep/external_data_initial_metadata_table.csv",
	output:
		filtered="data/external_sc_data_prep/sc_filtered.rds",
		meta="data/external_sc_data_prep/sc_metadata_filtered.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/filter_data.R"

rule write_final_metadata_table_ext_sc:
	input:
		meta="data/external_sc_data_prep/sc_metadata_filtered.csv",
	output:
		table="tables/external_sc_data_prep/external_data_final_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/write_metadata_table.R"

rule run_azimuth_ext_sc:
	input:
		data="data/external_sc_data_prep/sc_filtered.rds",
		table="tables/external_sc_data_prep/external_data_final_metadata_table.csv",
	output:
		annotated="data/external_sc_data_prep/sc_azimuth_annotated.rds",
		plot="plots/external_sc_data_prep/sc_azimuth_annotation.pdf",
	conda: "../env/azimuth.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/run_azimuth.R"

rule subset_samples_by_main_celltypes_ext_sc:
	input:
		data="data/external_sc_data_prep/sc_azimuth_annotated.rds",
	output:
		subset="data/external_sc_data_prep/sc_processed.rds",
		log="tables/external_sc_data_prep/sc_excluded_samples_log.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/subset_samples_by_main_celltypes.R"

rule summarize_final_data_included_ext_sc:
	input:
		sc="data/external_sc_data_prep/sc_processed.rds", 
	output:
		summary="tables/external_sc_data_prep/final_data_included_summary.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_sc_data_prep/summarize_final_data_included.R"
		

