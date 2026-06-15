rule run_deswan_immage:
	input:
		data="data/external_immage_data_prep/{mode}_pseudobulk_data_immage.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_immage.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan.R"

rule run_deswan_soundlife:
	input:
		data="data/external_soundlife_data_prep/{mode}_pseudobulk_data_soundlife.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_soundlife.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_soundlife.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_soundlife.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_soundlife.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan.R"

rule merge_immage_soundlife_combined1:
	input:
		"data/external_immage_data_prep/{mode}_pseudobulk_data_immage.csv",
		"data/external_soundlife_data_prep/{mode}_pseudobulk_data_soundlife.csv",
	output:
		"data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined1.csv",
		"tables/external_immage_analysis/{mode}_pseudobulk_data_column_summary_immage_soundlife_combined1.txt",
		"plots/external_immage_analysis/{mode}_pseudobulk_data_all_pca_immage_soundlife_combined1.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/merge_immage_soundlife_combined1.py"		

rule merge_immage_soundlife_combined2:
	input:
		"data/external_immage_data_prep/{mode}_pseudobulk_data_immage.csv",
		"data/external_soundlife_data_prep/{mode}_pseudobulk_data_soundlife.csv",
		"data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",	
	output:
		"data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined2.csv",
		"tables/external_immage_analysis/{mode}_pseudobulk_data_column_summary_immage_soundlife_combined2.txt",
		"plots/external_immage_analysis/{mode}_pseudobulk_data_all_pca_immage_soundlife_combined2.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/merge_immage_soundlife_combined2.py"		

rule merge_immage_soundlife_combined3:
	input:
		"data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",	
	output:
		"data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined3.csv",
		"tables/external_immage_analysis/{mode}_pseudobulk_data_column_summary_immage_soundlife_combined3.txt",
		"plots/external_immage_analysis/{mode}_pseudobulk_data_all_pca_immage_soundlife_combined3.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/merge_immage_soundlife_combined3.py"		

rule merge_immage_soundlife_combined4:
	input:
		"data/external_immage_data_prep/{mode}_pseudobulk_data_immage.csv",
		"data/external_soundlife_data_prep/{mode}_pseudobulk_data_soundlife.csv",
		"data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",	
	output:
		"data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined4.csv",
		"tables/external_immage_analysis/{mode}_pseudobulk_data_column_summary_immage_soundlife_combined4.txt",
		"plots/external_immage_analysis/{mode}_pseudobulk_data_all_pca_immage_soundlife_combined4.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/merge_immage_soundlife_combined4.py"		

rule run_deswan_immage_soundlife_combined1:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined1.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_immage_soundlife_combined1.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage_soundlife_combined1.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage_soundlife_combined1.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_soundlife_combined1.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan_combined.R"

rule run_deswan_immage_soundlife_combined2:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined2.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_immage_soundlife_combined2.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage_soundlife_combined2.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage_soundlife_combined2.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_soundlife_combined2.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan_combined.R"

rule run_deswan_immage_soundlife_combined3:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined3.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_immage_soundlife_combined3.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage_soundlife_combined3.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage_soundlife_combined3.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_soundlife_combined3.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan_combined.R"

rule run_deswan_immage_soundlife_combined4:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined4.csv",
	output:
		res="tables/external_immage_analysis/{mode}_deswan_res_immage_soundlife_combined4.rds",
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage_soundlife_combined4.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage_soundlife_combined4.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_soundlife_combined4.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan_combined.R"

rule run_deswan_immage_combined_loo:
	input:
		data="data/external_immage_data_prep/{mode}_pseudobulk_data_immage_combined.csv",
	output:
		coef="tables/external_immage_analysis/{mode}_deswan_coef_res_immage_combined_loo.csv",
		p="tables/external_immage_analysis/{mode}_deswan_p_res_immage_combined_loo.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_combined_loo.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/external_immage_analysis/run_deswan_combined_loo.R"

rule plot_deswan_immage:
	input:
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage.csv",
	output:
		plots="plots/external_immage_analysis/{mode}_deswan_q_res_by_gender_immage.pdf",
		plots2="plots/external_immage_analysis/{mode}_deswan_q_res_by_gender2_immage.pdf",
		plot="plots/external_immage_analysis/{mode}_deswan_q_res_all_immage.svg",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/plot_deswan.R"

rule plot_deswan_soundlife:
	input:
		q="tables/external_immage_analysis/{mode}_deswan_q_res_soundlife.csv",
	output:
		plots="plots/external_immage_analysis/{mode}_deswan_q_res_by_gender_soundlife.pdf",
		plots2="plots/external_immage_analysis/{mode}_deswan_q_res_by_gender2_soundlife.pdf",
		plot="plots/external_immage_analysis/{mode}_deswan_q_res_all_soundlife.svg",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/plot_deswan.R"

rule plot_deswan_immage_combined:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined{i}.csv",
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_soundlife_combined{i}.csv",
	output:
		plot1="plots/external_immage_analysis/{mode}_deswan_q_res_all_immage_soundlife_combined{i}.svg",
		plot2="plots/external_immage_analysis/{mode}_age_sex_distribution_immage_soundlife_combined{i}.svg",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/plot_deswan_immage_combined.R"

rule plot_deswan_immage_combined_loo:
	input:
		q="tables/external_immage_analysis/{mode}_deswan_q_res_immage_combined_loo.csv",
	output:
		plots="plots/external_immage_analysis/{mode}_deswan_q_res_combined_loo.pdf",
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/plot_deswan_combined_loo.R"

rule plot_sample_distribution_immage_soundlife_combined:
	input:
		data="data/external_immage_analysis/{mode}_pseudobulk_data_immage_soundlife_combined2.csv",
	output:
		plot="plots/external_immage_analysis/{mode}_immage_soundlife_combined2_sample_distribution.svg",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 40, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/external_immage_analysis/plot_sample_distribution.R"
