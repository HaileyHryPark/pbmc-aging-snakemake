rule get_full5ct_pseudobulk_mat_external_dis:
	input:
		data="data/external_dis_data_prep/{dataset}_{split}_processed.rds",
	output:
		pb="data/external_pseudobulk/{dataset}_{split}_full5ct_pseudobulk_data.csv",
	params: dataset="{dataset}"
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_pseudobulk/get_full5ct_pseudobulk_mat.R"

rule get_full5ct_pseudobulk_mat_external_sc:
	input:
		data="data/external_sc_data_prep/sc_processed.rds",
	output:
		pb="data/external_pseudobulk/sc_full5ct_pseudobulk_data.csv",
	params: dataset="sc"
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_pseudobulk/get_full5ct_pseudobulk_mat.R"

rule merge_full5ct_pseudobulk_data_external:
	input:
		"data/external_pseudobulk/sc_full5ct_pseudobulk_data.csv",
		expand("data/external_pseudobulk/ren_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 9)]), 
		expand("data/external_pseudobulk/wellcome_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 4)]), 
		expand("data/external_pseudobulk/combat_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 5)]), 
		"data/external_pseudobulk/ch_split01_full5ct_pseudobulk_data.csv",
		"data/external_pseudobulk/glaucoma_split01_full5ct_pseudobulk_data.csv",
		expand("data/external_pseudobulk/ra_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 2)]), 
		expand("data/external_pseudobulk/sle_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 6)]), 
	output:
		"data/external_pseudobulk/full5ct_pseudobulk_data_all.csv",
		"tables/external_pseudobulk/full5ct_pseudobulk_data_column_summary.txt",
		"plots/external_pseudobulk/full5ct_pseudobulk_data_all_pca.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "normal"
	script:
		"../scripts/external_pseudobulk/merge_full5ct_pseudobulk_data.py"

rule subset_allexp5ct_deswan_deg_pseudobulk_external:
	input:
		data="data/external_pseudobulk/full5ct_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/allexp5ct_deswan_q_deg_{gender}.csv"
	output:
		res="data/external_clock/allexp5ct_deswan_deg_pseudobulk_{gender}_data_all.csv",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/external_pseudobulk/subset_allexp5ct_deswan_deg_pseudobulk.R"

rule check_sample_distribution:
	input:
		data="data/external_pseudobulk/{mode}_pseudobulk_data_all.csv",
	output:
		age_res="tables/external_pseudobulk/{mode}_sample_distribution_age.txt",
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 40, walltime = "05:00:00", queue = "normal"
	script:
		"../scripts/external_pseudobulk/check_sample_distribution.R"

