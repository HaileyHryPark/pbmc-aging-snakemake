#rule run_deswan:
#	input:
#		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
#	output:
#		res="tables/internal_deswan/{mode}_deswan_res.rds",
#		coef="tables/internal_deswan/{mode}_deswan_coef_res.csv",
#		p="tables/internal_deswan/{mode}_deswan_p_res.csv",
#		q="tables/internal_deswan/{mode}_deswan_q_res.csv",
#	conda: "../env/internal_data_prep.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "normal"
#	script:
#		"../scripts/internal_deswan/run_deswan.R"
#
#rule run_deswan_by_dataset:
#	input:
#		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
#	output:
#		coef="tables/internal_deswan/{mode}_deswan_coef_res_by_dataset.csv",
#		p="tables/internal_deswan/{mode}_deswan_p_res_by_dataset.csv",
#		q="tables/internal_deswan/{mode}_deswan_q_res_by_dataset.csv",
#	conda: "../env/internal_data_prep.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "normal"
#	script:
#		"../scripts/internal_deswan/run_deswan_by_dataset.R"
#
#rule run_deswan_by_dataset_loo:
#	input:
#		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
#	output:
#		coef="tables/internal_deswan/{mode}_deswan_coef_res_by_dataset_loo.csv",
#		p="tables/internal_deswan/{mode}_deswan_p_res_by_dataset_loo.csv",
#		q="tables/internal_deswan/{mode}_deswan_q_res_by_dataset_loo.csv",
#	conda: "../env/internal_data_prep.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "normal"
#	script:
#		"../scripts/internal_deswan/run_deswan_by_dataset_loo.R"
#
#rule run_deswan_diff_params:
#	input:
#		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
#	output:
#		res="tables/internal_deswan/{mode}_deswan_q_res_by_diff_params.csv",
#	conda: "../env/internal_data_prep.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "normal"
#	script:
#		"../scripts/internal_deswan/run_deswan_diff_params.R"
#
#rule run_deswan_others:
#	input:
#		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
#	output:
#		res="tables/internal_deswan/{mode}_deswan_q_res_by_random_permutation.csv",
#	conda: "../env/internal_data_prep.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "normal"
#	script:
#		"../scripts/internal_deswan/run_deswan_others.R"
#
rule plot_deswan:
	input:
		q="tables/internal_deswan/{mode}_deswan_q_res.csv",
	output:
		plots="plots/internal_deswan/{mode}_deswan_q_res_by_gender.pdf",
		plots2="plots/internal_deswan/{mode}_deswan_q_res_by_gender2.pdf",
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan.R"

rule plot_deswan_by_dataset:
	input:
		q="tables/internal_deswan/{mode}_deswan_q_res_by_dataset.csv",
	output:
		plots="plots/internal_deswan/{mode}_deswan_q_res_by_gender_by_dataset.pdf",
	conda: "../env/internal_deswan.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan_diff_dataset.R"

rule plot_deswan_by_dataset_loo:
	input:
		q="tables/internal_deswan/{mode}_deswan_q_res_by_dataset_loo.csv",
	output:
		plots="plots/internal_deswan/{mode}_deswan_q_res_by_gender_by_dataset_loo.pdf",
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan_diff_dataset_loo.R"

rule plot_deswan_diff_params:
	input:
		res="tables/internal_deswan/{mode}_deswan_q_res_by_diff_params.csv",
	output:
		plot1="plots/internal_deswan/{mode}_deswan_q_res_by_diff_qvalues.pdf",
		plot2="plots/internal_deswan/{mode}_deswan_q_res_by_diff_buckets.pdf",
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan_diff_params.R"

rule plot_deswan_others:
	input:
		res="tables/internal_deswan/{mode}_deswan_q_res_by_random_permutation.csv",
	output:
		plot1="plots/internal_deswan/{mode}_deswan_q_res_by_random_permutation.pdf",
	conda: "../env/internal_deswan.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan_others.R"

rule get_deswan_deg:
	input:
		q="tables/internal_deswan/{mode}_deswan_q_res.csv",
	output:
		degb="tables/internal_deswan/{mode}_deswan_q_deg_both.csv",
		degf="tables/internal_deswan/{mode}_deswan_q_deg_female.csv",
		degm="tables/internal_deswan/{mode}_deswan_q_deg_male.csv",
	conda: "../env/internal_deswan.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/get_deswan_deg.R"	 

rule plot_deswan_deg_venn:
	input:
		degb="tables/internal_deswan/{mode}_deswan_q_deg_both.csv",
		degf="tables/internal_deswan/{mode}_deswan_q_deg_female.csv",
		degm="tables/internal_deswan/{mode}_deswan_q_deg_male.csv",
	output:
		venn="plots/internal_deswan/{mode}_deswan_q_deg_venn.pdf"
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_deswan/plot_deswan_deg_venn.R"	 

