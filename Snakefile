
configfile: "config.yaml"

rule all:
	input:
		## Internal data prep
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
		"tables/internal_data_prep/internal_data_final_metadata_table.csv",
		expand("data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt", dataset=config["internal_ds"]),
		"tables/internal_data_prep/final_data_included_summary.csv",
		"data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		"plots/internal_pseudobulk/full5ct_sample_distribution.pdf",
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv",

		## Internal DE-SWAN		
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender2.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset_loo.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_qvalues.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_buckets.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_random_permutation.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_deg_venn.pdf",

		## Internal clustering
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf", gender=["both","female","male"]),
		"plots/internal_clustering/allexp5ct_deswan_cluster_pie_by_celltype.pdf",
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_all_res.csv", gender=["both","female","male"]),
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_mitocarta_fa_res_plots.pdf",
		#"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.pdf",
		#"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_top_fa_score_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_all_res.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_loess.pdf",
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3.pdf", db=["GO","Reactome"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_{db}_network.pdf", db=["GO","Reactome"]),

		## Internal clock
		expand("plots/internal_clock/allexp5ct_deswan_deg_{model}_pred_comparison.pdf", model=["enet","xgboost","2dmlp"]),
		expand("plots/internal_clock/allexp5ct_deswan_deg_{model}_metric_comparison.pdf", model=["enet","xgboost","2dmlp"]),
		"plots/internal_clock/allexp5ct_deswan_deg_model_comparison_all.pdf",
		expand("tables/internal_clock/allexp5ct_deswan_deg_2dmlp_{gender}_model_shap_values.csv", gender=["both","female","male"]),
		
		## Internal celltype proportion analysis
		"data/internal_celltype_prop/celltype_prop_data_all.csv",
		"tables/internal_celltype_prop/celltype_prop_data_spearman_corr_age.csv",

		## External data analysis
		"tables/external_dis_data_prep/external_data_initial_metadata_table.csv",
		"tables/external_dis_data_prep/external_data_final_metadata_table.csv",
		expand("data/external_dis_data_prep/split_h5ad_by_sample_log_{dataset}.txt", dataset=["ren","wellcome","combat","ch","glaucoma","ra", "sle"]),
		"tables/external_dis_data_prep/final_data_included_summary.csv",
		"tables/external_sc_data_prep/final_data_included_summary.csv",
		"tables/external_pseudobulk/full5ct_sample_distribution_age.txt",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_metrics.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_scatter.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_disease.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_2dmlp_model_prediction_corrected_scatter.pdf",		

		### DNA methylation analysis
		#"data/dna_methylation/nsphs_norm.rds",		
		#"data/dna_methylation/nsphs_filtered.rds",
		#"tables/dna_methylation/nsphs_beta.csv",
		#"tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_both.rds",
		#"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_chromosome_pie.pdf",
		#"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_plots.pdf",
		#"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_line_vln_plots.pdf",

rule internal_prep:
	input:
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
		"tables/internal_data_prep/internal_data_final_metadata_table.csv",
		expand("data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt", dataset=config["internal_ds"]),
		"tables/internal_data_prep/final_data_included_summary.csv",
		"plots/internal_data_prep/final_data_cellcount_by_dataset.svg",
		"data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		"plots/internal_pseudobulk/full5ct_sample_distribution.pdf",
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv",

rule internal_deswan:
	input:
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender2.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset_loo.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_qvalues.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_buckets.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_random_permutation.pdf",
		"plots/internal_deswan/allexp5ct_deswan_q_deg_venn.pdf",

rule internal_clock:
	input:
		expand("plots/internal_clock/allexp5ct_deswan_deg_{model}_pred_comparison.pdf", model=["enet","xgboost","2dmlp"]),
		expand("plots/internal_clock/allexp5ct_deswan_deg_{model}_metric_comparison.pdf", model=["enet","xgboost","2dmlp"]),
		"plots/internal_clock/allexp5ct_deswan_deg_model_comparison_all.pdf",
		expand("tables/internal_clock/allexp5ct_deswan_deg_2dmlp_{gender}_model_shap_values.csv", gender=["both","female","male"]),
		"plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_pcaplots.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_corrected_shap_pcaplots.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_age_corrplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_age_corrplots_male.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_vlnplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_vlnplots_male.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_vlnplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_vlnplots_male.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_accel_vlnplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_accel_vlnplots_male.svg",
		
rule internal_clustering:
	input:
		#expand("tables/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf", gender=["both","female","male"]),
		"plots/internal_clustering/allexp5ct_deswan_cluster_pie_by_celltype.pdf",
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_all_res.csv", gender=["both","female","male"]),
		#"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_mitocarta_fa_res_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_all_res_network_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_all_res_network_plots2.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_top_fa_score_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_all_res.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_loess.pdf",
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3.pdf", db=["GO","Reactome"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_{db}_network.pdf", db=["GO","Reactome"]),

rule internal_clustering2:
	input:
		expand("plots/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_trajectories_{gender}.svg", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_multirun_cluster_trajectories_{gender}.svg", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}_final.svg", gender=["both","female","male"]),
		"plots/internal_clustering/allexp5ct_deswan_deg_gender_flowplot_by_merged_clusters.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_gender_venn_by_merged_clusters.pdf",	
		"plots/internal_clustering/allexp5ct_deswan_cluster_pie_by_celltype.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_GO_res_network3.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_GO_res_network3_fli.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_GO_res_network3_mef.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_loess.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_celltype_pie.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_all_res.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_GO_network_fcimei.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_female_specific_cluster_genes.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_fiu_cluster_genes.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_fth1.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_ftl.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_fcimei_cluster_genes.svg",
		"tables/internal_clustering/allexp5ct_deswan_deg_female_late_increase_ias_test.txt",

rule internal_cor:
	input:
		"tables/internal_correlation/allexp5ct_deswan_deg_limma_spearman_corr_gsea_all_res.csv",
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_spearman_corr_gsea_all_res_plots.pdf",
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_top_percentile_gender_corr_loess.pdf",
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_corr_comparison_specific_clusters.pdf",

rule internal_ctp:
	input:
		"data/internal_celltype_prop/celltype_prop_data_all.csv",
		"tables/internal_celltype_prop/celltype_prop_data_spearman_corr_age.csv",
		#"plots/internal_celltype_prop/celltype_prop_data_raw_scatter.pdf",
		#"plots/internal_celltype_prop/celltype_data_deswan_corr.pdf",

rule external:
	input:
		"tables/external_dis_data_prep/external_data_initial_metadata_table.csv",
		"tables/external_dis_data_prep/external_data_final_metadata_table.csv",
		expand("data/external_dis_data_prep/split_h5ad_by_sample_log_{dataset}.txt", dataset=["ren","wellcome","combat","ch","glaucoma","ra", "sle"]),
		"tables/external_dis_data_prep/final_data_included_summary.csv",
		"tables/external_sc_data_prep/final_data_included_summary.csv",
		"tables/external_pseudobulk/full5ct_sample_distribution_age.txt",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_metrics.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_scatter.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_disease.pdf",
		"plots/external_clock/allexp5ct_deswan_deg_2dmlp_model_prediction_corrected_scatter.pdf",

rule downstream:
	input:
		"plots/external_clock/allexp5ct_deswan_deg_2dmlp_model_prediction_corrected_scatter.pdf",		
		expand("plots/cluster_score/allexp5ct_deswan_deg_{cohort}_cluster_scores_{gender}_age_accel_radarplot.pdf", cohort=["internal","external"], gender=["female","male"]),

rule dnam:
	input:
		"data/dna_methylation/nsphs_norm.rds",		
		"data/dna_methylation/nsphs_filtered.rds",
		"tables/dna_methylation/nsphs_beta.csv",
		"plots/dna_methylation/nsphs_sample_distribution_age_gender.pdf",
		"tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_both.rds",
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_chromosome_pie.pdf",
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_plots.pdf",
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_line_vln_plots.pdf",
		expand("plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_line_vln_plot{i}.svg", i=["1","2","3","4","5"]),

rule cima_rna:
	input:
		"data/cima_rna_data_prep/split_h5ad_by_donor_log_cima_rna.txt",

		"tables/cima_rna_data_prep/final_data_included_summary.csv",
		"data/cima_rna_data_prep/allexpcima_pseudobulk_data_cima_rna.csv",
		"plots/cima_rna_data_prep/cima_rna_sample_distribution.svg",
		"plots/cima_rna_analysis/allexp5ct_deswan_q_res_by_gender.pdf",
		"plots/cima_rna_analysis/allexp5ct_deswan_q_deg_cima_rna_venn.pdf",
		"plots/cima_rna_analysis/allexp5ct_deswan_deg_model_comparison_metrics.pdf",
		"plots/cima_rna_analysis/allexp5ct_deswan_deg_2dmlp_model_prediction_corrected_scatter.pdf",

rule cima_atac:
	input:
		"data/cima_atac_data/cima_atac_metadata_prefilter.csv",

##### load rules #####

include: "rules/r_package_install.smk"
include: "rules/internal_data_prep.smk"
include: "rules/internal_pseudobulk.smk"
include: "rules/internal_deswan.smk"
include: "rules/internal_clustering.smk"
include: "rules/internal_correlation.smk"
include: "rules/internal_celltype_prop.smk"
include: "rules/internal_clock.smk"
include: "rules/external_dis_data_prep.smk"
include: "rules/external_sc_data_prep.smk"
include: "rules/external_pseudobulk.smk"
include: "rules/external_clock.smk"
include: "rules/cluster_score.smk"
include: "rules/dna_methylation.smk"
include: "rules/cima_rna_data_prep.smk"
include: "rules/cima_rna_analysis.smk"
include: "rules/cima_atac_data.smk"
