rule run_gender_corr_limma_deswan_deg_subset:
	input:
		limma="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}.csv",
	output:
		cor="tables/internal_correlation/{mode}_deswan_deg_limma_{gender}_corr.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "80:00:00", queue = "super"
	script:
		"../scripts/internal_correlation/run_gender_corr_limma_deswan_deg_subset.R"

rule run_mfuzz_merged_clusters_spearman_corr_gsea_all:
	input:
		both_cor="tables/internal_correlation/{mode}_deswan_deg_limma_both_corr.csv",
		female_cor="tables/internal_correlation/{mode}_deswan_deg_limma_female_corr.csv",
		male_cor="tables/internal_correlation/{mode}_deswan_deg_limma_male_corr.csv",
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		annot_cor="tables/internal_correlation/{mode}_deswan_deg_limma_corr_annotated_all.csv",
		res="tables/internal_correlation/{mode}_deswan_deg_limma_spearman_corr_gsea_all_res.rds",
		res_df="tables/internal_correlation/{mode}_deswan_deg_limma_spearman_corr_gsea_all_res.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "99:00:00", queue = "super"
	script:
		"../scripts/internal_correlation/run_mfuzz_merged_clusters_spearman_corr_gsea_all.R"

rule plot_top_percentile_gender_corr:
	input:
		annot_cor="tables/internal_correlation/{mode}_deswan_deg_limma_corr_annotated_all.csv",
		fit_res_f="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_female.csv",
		span_res_f="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_female.csv",
		fit_res_m="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_male.csv",
		span_res_m="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_male.csv",
	output:
		venn="plots/internal_correlation/{mode}_deswan_deg_limma_top_percentile_gender_corr_venn.pdf",
		pie="plots/internal_correlation/{mode}_deswan_deg_limma_top_percentile_gender_corr_pie.pdf",
		loess="plots/internal_correlation/{mode}_deswan_deg_limma_top_percentile_gender_corr_loess.pdf",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "80:00:00", queue = "super"
	script:
		"../scripts/internal_correlation/plot_top_percentile_gender_corr.R"

rule plot_mfuzz_merged_clusters_spearman_corr_gsea_all_network:
	input:
		res="tables/internal_correlation/{mode}_deswan_deg_limma_spearman_corr_gsea_all_res.csv",
	output:
		plots="plots/internal_correlation/{mode}_deswan_deg_limma_spearman_corr_gsea_all_res_plots.pdf",
	conda: "../env/internal_downstream2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_correlation/plot_mfuzz_merged_clusters_spearman_corr_gsea_all_network.R"

rule plot_corr_comparison_specific_clusters:
	input:
		annot_cor="tables/internal_correlation/{mode}_deswan_deg_limma_corr_annotated_all.csv",
	output:
		plots="plots/internal_correlation/{mode}_deswan_deg_limma_corr_comparison_specific_clusters.pdf",
	conda: "../env/internal_downstream2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_correlation/plot_corr_comparison_specific_clusters.R"

