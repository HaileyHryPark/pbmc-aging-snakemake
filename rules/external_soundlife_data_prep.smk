rule extract_metadata_ext_soundlife:
	input:
		"data/external_soundlife_data_prep/soundlife_{group}.h5ad",
	output:
		"data/external_soundlife_data_prep/soundlife_{group}_qced.h5ad",
		"data/external_soundlife_data_prep/soundlife_{group}_metadata.csv",
		"data/external_soundlife_data_prep/soundlife_{group}_ensembl_to_symbol.csv",
	params: dataset="soundlife", group="{group}"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/extract_metadata.py"

rule write_initial_metadata_table_ext_soundlife:
	input:
		ofcn="data/external_soundlife_data_prep/soundlife_OldF_CMVn_metadata.csv",
		ofcp="data/external_soundlife_data_prep/soundlife_OldF_CMVp_metadata.csv",
		omcn="data/external_soundlife_data_prep/soundlife_OldM_CMVn_metadata.csv",
		omcp="data/external_soundlife_data_prep/soundlife_OldM_CMVp_metadata.csv",
		yfcn="data/external_soundlife_data_prep/soundlife_YoungF_CMVn_metadata.csv",
		yfcp="data/external_soundlife_data_prep/soundlife_YoungF_CMVp_metadata.csv",
		ymcn="data/external_soundlife_data_prep/soundlife_YoungM_CMVn_metadata.csv",
		ymcp="data/external_soundlife_data_prep/soundlife_YoungM_CMVp_metadata.csv",
	output:
		table="tables/external_soundlife_data_prep/external_data_initial_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/write_metadata_table.R"

rule filter_data_ext_soundlife:
	input:
		"data/external_soundlife_data_prep/soundlife_{group}_qced.h5ad",
		"tables/external_soundlife_data_prep/external_data_initial_metadata_table.csv",
	output:
		"data/external_soundlife_data_prep/soundlife_{group}_filtered.h5ad",
		"data/external_soundlife_data_prep/soundlife_{group}_metadata_filtered.csv",
	params: dataset="soundlife"
	conda: "../env/internal_data_prep_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 300, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/filter_data.py"

rule write_final_metadata_table_ext_soundlife:
	input:
		ofcn="data/external_soundlife_data_prep/soundlife_OldF_CMVn_metadata_filtered.csv",
		ofcp="data/external_soundlife_data_prep/soundlife_OldF_CMVp_metadata_filtered.csv",
		omcn="data/external_soundlife_data_prep/soundlife_OldM_CMVn_metadata_filtered.csv",
		omcp="data/external_soundlife_data_prep/soundlife_OldM_CMVp_metadata_filtered.csv",
		yfcn="data/external_soundlife_data_prep/soundlife_YoungF_CMVn_metadata_filtered.csv",
		yfcp="data/external_soundlife_data_prep/soundlife_YoungF_CMVp_metadata_filtered.csv",
		ymcn="data/external_soundlife_data_prep/soundlife_YoungM_CMVn_metadata_filtered.csv",
		ymcp="data/external_soundlife_data_prep/soundlife_YoungM_CMVp_metadata_filtered.csv",
	output:
		table="tables/external_soundlife_data_prep/external_data_final_metadata_table.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 380, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/write_metadata_table.R"

# manual_seu_conversion.R
rule convert_h5ad_to_seu_ext_soundlife:
	input:
		data="data/external_soundlife_data_prep/soundlife_{group}_filtered.h5ad",
	output:
		seu="data/external_soundlife_data_prep/soundlife_{group}_filtered.rds",
	conda: "../env/anndata.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/convert_h5ad_to_seu.R"

## Had to install manually and use server singularity (server outdated - old glibc); --use-singularity
rule run_azimuth_by_split_ext_soundlife:
	input:
		data="data/external_soundlife_data_prep/soundlife_{group}_filtered.rds",
	output:
		annotated="data/external_soundlife_data_prep/soundlife_{group}_azimuth_annotated.rds",
		plot="plots/external_soundlife_data_prep/soundlife_{group}_azimuth_annotation.pdf",
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/run_azimuth_by_split.R"

rule subset_samples_by_main_celltypes_ext_soundlife:
	input:
		data="data/external_soundlife_data_prep/soundlife_{group}_azimuth_annotated.rds",
	output:
		subset="data/external_soundlife_data_prep/soundlife_{group}_processed.rds",
		log="tables/external_soundlife_data_prep/soundlife_{group}_excluded_samples_log.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/subset_samples_by_main_celltypes.R"

rule summarize_final_data_included_ext_soundlife:
	input:
		data=expand("data/external_soundlife_data_prep/soundlife_{group}_processed.rds", group=["OldF_CMVn", "OldF_CMVp", "OldM_CMVn", "OldM_CMVp", "YoungF_CMVn", "YoungF_CMVp", "YoungM_CMVn", "YoungM_CMVp"]),
	output:
		summary="tables/external_soundlife_data_prep/final_data_included_summary.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 200, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/summarize_final_data_included.R"
		
rule get_full5ct_pseudobulk_mat_soundlife:
	input:
		data="data/external_soundlife_data_prep/soundlife_{group}_processed.rds",
	output:
		pb="data/external_soundlife_data_prep/soundlife_{group}_full5ct_pseudobulk_data.csv",
	params: dataset="soundlife"
	conda: "../env/internal_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/get_full5ct_pseudobulk_mat.R"

rule merge_full5ct_pseudobulk_data_soundlife:
	input:
		expand("data/external_soundlife_data_prep/soundlife_{group}_full5ct_pseudobulk_data.csv", group=["OldF_CMVn", "OldF_CMVp", "OldM_CMVn", "OldM_CMVp", "YoungF_CMVn", "YoungF_CMVp", "YoungM_CMVn", "YoungM_CMVp"]),
	output:
		"data/external_soundlife_data_prep/full5ct_pseudobulk_data_soundlife.csv",
		"tables/external_soundlife_data_prep/full5ct_pseudobulk_data_soundlife_column_summary.txt",
		"plots/external_soundlife_data_prep/full5ct_pseudobulk_data_soundlife_pca.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/merge_full5ct_pseudobulk_data.py"

rule subset_allexp5ct_pseudobulk_soundlife:
	input:
		"data/external_soundlife_data_prep/full5ct_pseudobulk_data_soundlife.csv",
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv"
	output:
		"data/external_soundlife_data_prep/allexp5ct_pseudobulk_data_soundlife.csv",
		"data/external_soundlife_data_prep/allexp5ct_pseudobulk_data_soundlife_combined.csv",
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_soundlife_data_prep/subset_allexp5ct_pseudobulk.py"

rule plot_sample_distribution_soundlife:
        input:
                data="data/external_soundlife_data_prep/full5ct_pseudobulk_data_soundlife.csv",
        output:
                plot="plots/external_soundlife_data_prep/soundlife_sample_distribution.svg",
        conda: "../env/internal_downstream.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 40, walltime = "05:00:00", queue = "super"
        script:
                "../scripts/external_soundlife_data_prep/plot_sample_distribution.R"

