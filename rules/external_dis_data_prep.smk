rule extract_metadata_ext_dis:
	input:
		"data/external_dis_data_prep/{dataset}.h5ad",
	output:
		"data/external_dis_data_prep/{dataset}_qced.h5ad",
		"data/external_dis_data_prep/{dataset}_metadata.csv",
		"data/external_dis_data_prep/{dataset}_ensembl_to_symbol.csv",
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/extract_metadata.py"

rule write_initial_metadata_table_ext_dis:
	input:
		ren="data/external_dis_data_prep/ren_metadata.csv",
		wellcome="data/external_dis_data_prep/wellcome_metadata.csv",
		combat="data/external_dis_data_prep/combat_metadata.csv",
		ch="data/external_dis_data_prep/ch_metadata.csv",
		glaucoma="data/external_dis_data_prep/glaucoma_metadata.csv",
		ra="data/external_dis_data_prep/ra_metadata.csv",
	output:
		table="tables/external_dis_data_prep/external_data_initial_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/write_metadata_table.R"

rule filter_data_ext_dis:
	input:
		"data/external_dis_data_prep/{dataset}_qced.h5ad",
		"tables/external_dis_data_prep/external_data_initial_metadata_table.csv",
	output:
		"data/external_dis_data_prep/{dataset}_filtered.h5ad",
		"data/external_dis_data_prep/{dataset}_metadata_filtered.csv",
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/filter_data.py"

rule write_final_metadata_table_ext_dis:
	input:
		ren="data/external_dis_data_prep/ren_metadata_filtered.csv",
		wellcome="data/external_dis_data_prep/wellcome_metadata_filtered.csv",
		combat="data/external_dis_data_prep/combat_metadata_filtered.csv",
		ch="data/external_dis_data_prep/ch_metadata_filtered.csv",
		glaucoma="data/external_dis_data_prep/glaucoma_metadata_filtered.csv",
		ra="data/external_dis_data_prep/ra_metadata_filtered.csv",
	output:
		table="tables/external_dis_data_prep/external_data_final_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/write_metadata_table.R"

rule split_h5ad_by_sample_ext_dis:
	input:
		"data/external_dis_data_prep/{dataset}_filtered.h5ad",
	output:
		"data/external_dis_data_prep/split_h5ad_by_sample_log_{dataset}.txt"
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/split_h5ad_by_sample.py"

rule split_h5ad_by_sample_sle:
	input:
		"data/internal_data_prep/perez_filtered_others.h5ad",
	output:
		"data/external_dis_data_prep/split_h5ad_by_sample_log_sle.txt"
	params: dataset="sle"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/split_h5ad_by_sample_sle.py"

# manual_seu_conversion.R

rule run_azimuth_by_split_ext_dis:
	input:
		data="data/external_dis_data_prep/{dataset}_filtered_{split}.rds",
	output:
		annotated="data/external_dis_data_prep/{dataset}_{split}_azimuth_annotated.rds",
		plot="plots/external_dis_data_prep/{dataset}_{split}_azimuth_annotation.pdf",
	conda: "../env/azimuth.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/run_azimuth_by_split.R"

rule subset_samples_by_main_celltypes_ext_dis:
	input:
		data="data/external_dis_data_prep/{dataset}_{split}_azimuth_annotated.rds",
	output:
		subset="data/external_dis_data_prep/{dataset}_{split}_processed.rds",
		log="tables/external_dis_data_prep/{dataset}_{split}_excluded_samples_log.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/subset_samples_by_main_celltypes.R"

rule summarize_final_data_included_ext_dis:
	input:
		ren=expand("data/external_dis_data_prep/ren_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 9)]),
		wellcome=expand("data/external_dis_data_prep/wellcome_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 4)]),
		combat=expand("data/external_dis_data_prep/combat_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 5)]),
		ch="data/external_dis_data_prep/ch_split01_processed.rds", 
		glaucoma="data/external_dis_data_prep/glaucoma_split01_processed.rds",
		ra=expand("data/external_dis_data_prep/ra_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 2)]),
		sle=expand("data/external_dis_data_prep/sle_{split}_processed.rds", split=[f"split{i:02d}" for i in range(1, 6)]),
	output:
		summary="tables/external_dis_data_prep/final_data_included_summary.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_dis_data_prep/summarize_final_data_included.R"
		

