rule run_limma_deswan_pseudobulk_data:
	input:
		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/{mode}_deswan_q_deg_{gender}.csv",
	output:
		res="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}.csv",
	params: gender="{gender}"
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_limma_removebatcheffect_pseudobulk_data.R"

rule fit_gender_loess_zscore_limma_deswan_deg_subset:
	input:
		limma="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}.csv",
	output:
		zscaled="tables/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}_zscaled.csv",
		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_{gender}.csv",
		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_{gender}.csv",
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "80:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/fit_gender_loess_zscore_limma_deswan_deg_subset.R"

rule estimate_params_mfuzz_loess_fitted:
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.pdf",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/estimate_params_mfuzz_loess_fitted.R"

rule cluster_mfuzz_loess_fitted:
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
	output:
		cnumplot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_elbowplot_{gender}.pdf",
		corrplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_initial_{gender}.pdf",
		mfuzzplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_initial_{gender}.pdf",
		corrplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_merged_{gender}.pdf",
		mfuzzplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_merged_{gender}.pdf",
		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/cluster_mfuzz_loess_fitted.R"

rule plot_mfuzz_merged_clusters_loess:
	input:
		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_{gender}.csv",
		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_{gender}.csv",
		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv",
	output:
		annotated="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf",
	params: gender="{gender}"
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_loess.R"

rule plot_mfuzz_merged_clusters_info:
	input:
		both="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		flowplot="plots/internal_clustering/{mode}_deswan_deg_gender_flowplot_by_merged_clusters.pdf",
		venn="plots/internal_clustering/{mode}_deswan_deg_gender_venn_by_merged_clusters.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_info.R"

rule run_mfuzz_merged_clusters_fa:
	input:
		df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_{gender}_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_{gender}_res_no_ribo.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_fa.R"

rule run_mfuzz_merged_clusters_msigdb_fa:
	input:
		df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_{gender}_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_{gender}_res_no_ribo.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_msigdb_fa.R"

rule run_mfuzz_merged_clusters_disgenet2r:
	input:
		df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res_no_ribo.csv",
		plot1="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_plot.pdf",
		plot2="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_plot_no_ribo.pdf",
	conda: "../env/azimuth.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_disgenet2r.R"

rule plot_mfuzz_merged_clusters_fa:
	input:
		both1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_both_res.csv",
		both2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_both_res_no_ribo.csv",
		female1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_female_res.csv",
		female2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_female_res_no_ribo.csv",
		male1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_male_res.csv",
		male2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_male_res_no_ribo.csv",
	output:
		plot1="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_plot.pdf",
		plot2="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_plot_no_ribo.pdf",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa.R"

