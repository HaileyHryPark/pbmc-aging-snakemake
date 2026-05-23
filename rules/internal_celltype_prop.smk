rule get_celltype_prop_mat:
	input:
		data="data/internal_data_prep/{dataset}_{split}_processed.rds",
	output:
		ctp="data/internal_celltype_prop/{dataset}_{split}_celltype_prop_data.csv",
		meta="data/internal_celltype_prop/{dataset}_{split}_celltype_meta_data.csv",
		meta2="data/internal_celltype_prop/{dataset}_{split}_celltype_meta_data_summary.csv",
	params: dataset="{dataset}"
	conda: "../env/external_dis_data_prep.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/get_celltype_prop_mat.R"

rule merge_celltype_prop_data:
	input:
		expand("data/internal_celltype_prop/onek1k_{split}_celltype_prop_data.csv", split=[f"split{i:02d}" for i in range(1, 11)]), 
		expand("data/internal_celltype_prop/aida_{split}_celltype_prop_data.csv", split=[f"split{i:02d}" for i in range(1, 22)]), 
		expand("data/internal_celltype_prop/perez_{split}_celltype_prop_data.csv", split=[f"split{i:02d}" for i in range(1, 5)]), 
		expand("data/internal_celltype_prop/marina_{split}_celltype_prop_data.csv", split=[f"split{i:02d}" for i in range(1, 7)]), 
	output:
		"data/internal_celltype_prop/celltype_prop_data_all.csv",
		"tables/internal_celltype_prop/celltype_prop_data_column_summary.txt",
		"plots/internal_celltype_prop/celltype_prop_data_all_pca.pdf"
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/merge_celltype_prop_data.py"

rule run_limma_celltype_prop_data:
	input:
		data="data/internal_celltype_prop/celltype_prop_data_all.csv",
	output:
		res="data/internal_celltype_prop/celltype_prop_data_limma.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/run_limma_removebatcheffect_celltype_prop_data.R"

rule plot_celltype_prop_scatter:
	input:
		raw="data/internal_celltype_prop/celltype_prop_data_all.csv",
		data="data/internal_celltype_prop/celltype_prop_data_limma.csv",
	output:
		table="tables/internal_celltype_prop/celltype_prop_data_spearman_corr_age.csv",
		plot1="plots/internal_celltype_prop/celltype_prop_data_limma_scatter.pdf",
		plot2="plots/internal_celltype_prop/celltype_prop_data_raw_scatter.pdf",
		plot3="plots/internal_celltype_prop/celltype_prop_data_raw_violin.pdf",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/plot_celltype_prop_data.R"

rule merge_celltype_meta_data:
	input:
		expand("data/internal_celltype_prop/onek1k_{split}_celltype_meta_data.csv", split=[f"split{i:02d}" for i in range(1, 11)]), 
		expand("data/internal_celltype_prop/aida_{split}_celltype_meta_data.csv", split=[f"split{i:02d}" for i in range(1, 22)]), 
		expand("data/internal_celltype_prop/perez_{split}_celltype_meta_data.csv", split=[f"split{i:02d}" for i in range(1, 5)]), 
		expand("data/internal_celltype_prop/marina_{split}_celltype_meta_data.csv", split=[f"split{i:02d}" for i in range(1, 7)]), 
	output:
		"data/internal_celltype_prop/celltype_meta_data_all.csv",
		"tables/internal_celltype_prop/celltype_meta_data_column_summary.txt",
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/merge_celltype_meta_data.py"

rule plot_celltype_data_deswan_corr:
	input:
		degb="tables/internal_deswan/allexp5ct_deswan_q_deg_both.csv",
		degf="tables/internal_deswan/allexp5ct_deswan_q_deg_female.csv",
		degm="tables/internal_deswan/allexp5ct_deswan_q_deg_male.csv",
		table="data/internal_celltype_prop/celltype_meta_data_all.csv",
	output:
		plot="plots/internal_celltype_prop/celltype_data_deswan_corr.pdf",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/plot_celltype_data_deswan_corr.R"

rule run_celltype_prop_propeller:
	input:
		data="data/internal_celltype_prop/celltype_meta_data_all.csv",
	output:
		l1res_anova="tables/internal_celltype_prop/l1_propeller_{gender}_anova.csv",
		l1res_age="tables/internal_celltype_prop/l1_propeller_{gender}_age.csv",
		l1_vlnplot="plots/internal_celltype_prop/l1_propeller_{gender}_vlnplot.svg",
		l2res_anova="tables/internal_celltype_prop/l2_propeller_{gender}_anova.csv",
		l2res_age="tables/internal_celltype_prop/l2_propeller_{gender}_age.csv",
		l2_vlnplot="plots/internal_celltype_prop/l2_propeller_{gender}_vlnplot.svg",
	params: gender="{gender}"
	conda: "../env/internal_celltype_prop.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_celltype_prop/run_celltype_prop_propeller.R"

