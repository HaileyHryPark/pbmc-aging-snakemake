rule run_limma_deswan_pseudobulk_data:
	input:
		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
		degb="tables/internal_deswan/{mode}_deswan_q_deg_both.csv",
                degf="tables/internal_deswan/{mode}_deswan_q_deg_female.csv",
                degm="tables/internal_deswan/{mode}_deswan_q_deg_male.csv",
	output:
		res="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_all.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/run_limma_removebatcheffect_pseudobulk_data.R"

rule fit_gender_loess_zscore_limma_deswan_deg_subset:
	input:
		limma="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_all.csv",
	output:
		zscaled="tables/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_all_zscaled.csv",
		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted.csv",
		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res.csv",
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat.txt",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "80:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/fit_gender_loess_zscore_limma_deswan_deg_subset.R"

rule estimate_params_mfuzz_loess_fitted:
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat.txt",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num.pdf",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/estimate_params_mfuzz_loess_fitted.R"

rule cluster_mfuzz_loess_fitted:
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat.txt",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num.csv",
	output:
		cnumplot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_elbowplot.pdf",
		corrplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_initial.pdf",
		mfuzzplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_initial.pdf",
		corrplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_merged.pdf",
		mfuzzplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_merged.pdf",
		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/internal_clustering/cluster_mfuzz_loess_fitted.R"

#rule plot_mfuzz_clusters_heatmap:
#	input:
#		fit="tables/internal_clustering/{mode}_deswan_deg_loess_fitted.csv",
#		deg="tables/internal_clustering/{mode}_pseudobulk_data_deswan_q_deg.csv",
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
#	output:
#		plots="plots/internal_clustering/{mode}_deswan_deg_mfuzz_cluster_celltype_loess_heatmap.pdf",
#	conda: "../env/internal_clustering.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_clusters_heatmap.R"
#
#rule plot_mfuzz_initial_clusters_heatmap:
#	input:
#		fit="tables/internal_clustering/{mode}_deswan_deg_loess_fitted.csv",
#		deg="tables/internal_clustering/{mode}_pseudobulk_data_deswan_q_deg.csv",
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
#	output:
#		plot1="plots/internal_clustering/{mode}_deswan_deg_mfuzz_initial_cluster_celltype_loess_heatmap_all_ct.pdf",
#		plot2="plots/internal_clustering/{mode}_deswan_deg_mfuzz_initial_cluster_celltype_loess_heatmap.pdf",
#	params:
#		mode="{mode}"
#	conda: "../env/internal_clustering.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_initial_clusters_heatmap.R"
#
#rule plot_mfuzz_clusters_pca:
#	input:
#		fit="tables/internal_clustering/{mode}_deswan_deg_loess_fitted.csv",
#		deg="tables/internal_clustering/{mode}_pseudobulk_data_deswan_q_deg.csv",
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
#	output:
#		var_clust_df_annotated="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_annotated.csv",
#		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_clusters_pca.pdf",
#	params:
#		mode="{mode}"
#	conda: "../env/combined_pseudobulk_deg.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_clusters_pca.R"
#
#rule plot_mfuzz_clusters_umap:
#	input:
#		"tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_annotated.csv",
#		"tables/internal_clustering/{mode}_deswan_deg_loess_fitted.csv",
#	output:
#		"plots/internal_clustering/{mode}_deswan_deg_mfuzz_clusters_umap.pdf",
#	conda: "../env/combined_pseudobulk_py.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_clusters_umap.py"
#
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
#rule plot_mfuzz_initial_clusters_info:
#	input:
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment.csv",
#		deg="tables/internal_clustering/{mode}_pseudobulk_data_deswan_q_deg.csv",
#	output:
#		flowplot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_initial_cluster_gender_flowplot_by_celltype.pdf",
#		barplot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_initial_cluster_celltype_barplot_by_cluster.pdf",
#		piechart="plots/internal_clustering/{mode}_deswan_deg_mfuzz_initial_cluster_piechart_by_cluster.pdf",
#	conda: "../env/internal_clustering.yaml"
#	params:
#		mode="{mode}"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "normal"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_initial_clusters_info.R"
#
