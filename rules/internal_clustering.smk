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

#rule plot_mfuzz_clusters_info:
#	input:
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
#		deg="tables/internal_clustering/{mode}_pseudobulk_data_deswan_q_deg.csv",
#	output:
#		flowplot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_cluster_gender_flowplot_by_celltype.pdf",
#		barplot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_cluster_celltype_barplot_by_cluster.pdf",
#		piechart="plots/internal_clustering/{mode}_deswan_deg_mfuzz_cluster_piechart_by_cluster.pdf",
#	conda: "../env/internal_clustering.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_clusters_info.R"
#

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

rule run_mfuzz_merged_clusters_disgenet2r:
	input:
		df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res_no_ribo.csv",
		plot1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_plot.pdf",
		plot2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_plot_no_ribo.pdf",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_disgenet2r.R"
