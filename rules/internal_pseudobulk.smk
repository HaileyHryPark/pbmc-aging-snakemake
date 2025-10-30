rule get_full5ct_pseudobulk_mat:
       input:
               data="data/internal_data_prep/{dataset}_{split}_processed.rds",
       output:
               pb="data/internal_pseudobulk/{dataset}_{split}_full5ct_pseudobulk_data.csv",
       params: dataset="{dataset}"
       conda: "../env/internal_data_prep.yaml"
       threads: 1
       resources: ngpus = 0, mem_gb = 100, walltime = "20:00:00", queue = "normal"
       script:
               "../scripts/internal_pseudobulk/get_full5ct_pseudobulk_mat.R"

rule merge_full5ct_pseudobulk_data:
       input:
               expand("data/internal_pseudobulk/onek1k_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 10)]), 
               expand("data/internal_pseudobulk/aida_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 21)]), 
               expand("data/internal_pseudobulk/perez_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 4)]), 
               expand("data/internal_pseudobulk/marina_{split}_full5ct_pseudobulk_data.csv", split=[f"split{i:02d}" for i in range(1, 6)]), 
       output:
               "data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
               "tables/internal_pseudobulk/full5ct_pseudobulk_data_column_summary.txt",
               "plots/internal_pseudobulk/full5ct_pseudobulk_data_all_pca.pdf"
       conda: "../env/internal_pseudobulk_py.yaml"
       threads: 1
       resources: ngpus = 0, mem_gb = 100, walltime = "10:00:00", queue = "normal"
       script:
               "../scripts/internal_pseudobulk/merge_full5ct_pseudobulk_data.py"

rule subset_allexp5ct_pseudobulk:
	input:
		"data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
	output:
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv",
	conda: "../env/internal_pseudobulk_py.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "normal"
	script:
		"../scripts/internal_pseudobulk/subset_allexp5ct_pseudobulk.py"

rule plot_sample_distribution:
        input:
                data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
        output:
                plot1="plots/internal_pseudobulk/{mode}_sample_distribution.pdf",
                plot2="plots/internal_pseudobulk/{mode}_sample_distribution2.pdf",
                plot_ds="plots/internal_pseudobulk/{mode}_sample_distribution_by_dataset.pdf",
	conda: "../env/external_dis_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 40, walltime = "05:00:00", queue = "normal"
        script:
                "../scripts/internal_pseudobulk/plot_sample_distribution.R"

