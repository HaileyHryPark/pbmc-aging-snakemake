
configfile: "config.yaml"

rule all:
	input:
		## Internal data prep
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
		"tables/internal_data_prep/internal_data_final_metadata_table.csv",
		expand("data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt", dataset=config["internal_ds"]),
		"tables/internal_data_prep/final_data_included_summary.csv",
		"data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		# Fig. 1c
		"plots/internal_data_prep/all_umap.svg",
		# Supp Fig. 1
		"plots/internal_data_prep/all_umap_split_by_dataset.svg",
		# Fig. 1b
		"plots/internal_pseudobulk/full5ct_sample_distribution.pdf",
		# Extended Data Fig. 1b
		"plots/internal_pseudobulk/full5ct_sample_distribution_by_dataset.pdf",
		# Extended Data Fig. 1a
		"plots/internal_pseudobulk/full5ct_sample_distribution2.pdf",
		"data/internal_pseudobulk/allexp5ct_pseudobulk_data_all.csv",
		# Supp Fig. 4
		"plots/external_immage_data_prep/immage_sample_distribution.svg",
		"plots/external_soundlife_data_prep/soundlife_sample_distribution.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined1.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined2.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined3.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined4.svg",
		"plots/external_immage_analysis/allexp5ct_immage_soundlife_combined2_sample_distribution.svg",	
		# Supp Fig. 8
		"plots/internal_pseudobulk/cellnumber_pc1_correlation.svg",
		"plots/internal_pseudobulk/cellnumber_totalexp_correlation.svg",

		## Internal celltype proportion analysis
		# Table 1
		expand("tables/internal_celltype_prop/l1_propeller_{gender}_age.csv", gender=["both","female","male"]),
		expand("tables/internal_celltype_prop/l2_propeller_{gender}_age.csv", gender=["both","female","male"]),
		# Supp Fig. 2-3
		"plots/internal_celltype_prop/celltype_prop_data_raw_scatter.pdf",

		## Internal DE-SWAN		
		# Fig. 1e
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender.pdf",
		# Fig. 1f, Extended Data Fig. 2e
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender2.pdf",
		# Extended Data Fig. 1c
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset.pdf",
		# Extended Data Fig. 1d
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender_by_dataset_loo.pdf",
		# Extended Data Fig. 2a
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_qvalues.pdf",
		# Extended Data Fig. 2b
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_diff_buckets.pdf",
		# Extended Data Fig. 2c
		"plots/internal_deswan/allexp5ct_deswan_res_downsample.svg",
		# Extended Data Fig. 2d
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_random_permutation.pdf",
		# Fig. 1d
		"plots/internal_deswan/allexp5ct_deswan_q_deg_venn.pdf",

		## Internal clustering
		# Supp Fig. 5
		expand("plots/internal_clustering/allexp5ct_deswan_deg_pseudobulk_data_limma_{gender}_pca.svg", gender=["both","female","male"]),
		# Supp Fig. 6
		expand("plots/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_trajectories_{gender}.svg", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_multirun_cluster_trajectories_{gender}.svg", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_multirun_cluster_num_{gender}.svg", gender=["both","female","male"]),
		# Fig. 2a
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}_final.svg", gender=["both","female","male"]),
		# Extended Data Fig. 3a 
		"plots/internal_clustering/allexp5ct_deswan_cluster_pie_by_celltype.pdf",
		# Extended Data Fig. 3b
		"plots/internal_clustering/allexp5ct_deswan_deg_gender_flowplot_by_merged_clusters.pdf",
		# Extended Data Fig. 3c
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_top_percentile_gender_corr_bar.svg",
		# Fig. 2b
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_spearman_corr_gsea_all_res_plots.pdf",
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv", gender=["both","female","male"]),
		# Extended Data Fig. 4a
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.svg",
		# Fig. 3b,c
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_loess.pdf",
		# Fig. 3d-f
		"plots/internal_correlation/allexp5ct_deswan_deg_limma_corr_comparison_specific_clusters.pdf",
		# Extended Data Fig. 5a,c
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_GO_res_network3.pdf",
		# Fig. 3g,i, Extended Data Fig. 5f,h
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_GO_res_network3_{cl}.pdf", cl=["fci","fli","fiu","mef"]),
		# Fig. 3j
		expand("plots/internal_clustering/allexp5ct_deswan_deg_{gene}.svg", gene=["fth1","ftl"]),
		# Fig. 3l
		"plots/internal_clustering/allexp5ct_deswan_deg_fiu_cluster_genes.svg",
		# Extended Data Fig. 5b,e,g, Fig. 4n
		"plots/internal_clustering/allexp5ct_deswan_deg_female_specific_cluster_genes.svg",
		# Extended Data Fig. 5d
		"plots/internal_clustering/nhanes_serum_ferritin_plot.svg",
		# Fig. 4b
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.svg",

		### DNA methylation analysis
		"data/dna_methylation/nsphs_norm.rds",		
		"data/dna_methylation/nsphs_filtered.rds",
		"tables/dna_methylation/nsphs_beta.csv",
		"tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_both.rds",
		# Extended Data Fig. 6a
		"plots/dna_methylation/nsphs_sample_distribution_age_gender.pdf",
		# Fig. 4g-h
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_chromosome_pie.pdf",
		# Fig. 4d, i, Extended Data Fig. 6b
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_plots.pdf",
		# Fig. 4e-f, j-m, Extended Data Fig. 6c-f
		"plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_line_vln_plots.pdf",

		## Aging clock
		# Extended Data Fig. 7a-b, e-f
		"plots/internal_clock/allexp5ct_deswan_deg_model_comparison_all.pdf",
		# Supp Fig. 7
		"plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_pcaplots.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_corrected_shap_pcaplots.svg",
		# Fig. 5d-e
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_vlnplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_ag_vlnplots_male.svg",
		# Extended Data Fig. 8c-d
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_vlnplots_female.svg",
                "plots/internal_clock/allexp5ct_deswan_deg_2dmlp_shap_clust_vlnplots_male.svg",
		
		## External data analysis
		"data/external_pseudobulk/full5ct_pseudobulk_data_all.csv",
		# Fig. 5b-c, Extended Data Fig. 8a-b
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_scatter.pdf",
		# Extended Data Fig. 7c-d, g-h
		"plots/external_clock/allexp5ct_deswan_deg_model_comparison_metrics.pdf",
		# Extended Data Fig. 7i
		"plots/external_clock/allexp5ct_deswan_deg_mlp_comparison_healthy_disease.svg",

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
		#"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.svg",
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
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.svg",
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
		"plots/internal_clustering/nhanes_serum_ferritin_plot.svg",
		"tables/internal_clustering/allexp5ct_deswan_deg_female_late_increase_geneset_test.txt",
		"tables/internal_clustering/allexp5ct_deswan_deg_female_inverted_ushape_geneset_test.txt",
		"tables/internal_clustering/allexp5ct_deswan_deg_sex_specific_clusters_test.txt",

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
		"plots/internal_celltype_prop/celltype_prop_data_raw_violin.pdf",
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
		"plots/external_clock/allexp5ct_deswan_deg_mlp_comparison_healthy_disease.svg",
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

rule immage:
	input:
		"data/external_immage_data_prep/split_h5ad_by_sample_log_immage.txt",
		"tables/external_immage_data_prep/final_data_included_summary.csv",
		"data/external_immage_data_prep/allexp5ct_pseudobulk_data_immage.csv",
		"data/external_immage_data_prep/allexp5ct_pseudobulk_data_immage_combined.csv",
		"plots/external_immage_data_prep/immage_sample_distribution.svg",
		"plots/external_soundlife_data_prep/soundlife_sample_distribution.svg",
		"data/external_immage_analysis/allexp5ct_pseudobulk_data_immage_soundlife_combined1.csv",	
		"data/external_immage_analysis/allexp5ct_pseudobulk_data_immage_soundlife_combined2.csv",	
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_soundlife.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined1.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined2.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined3.svg",
		"plots/external_immage_analysis/allexp5ct_deswan_q_res_all_immage_soundlife_combined4.svg",
		"plots/external_immage_analysis/allexp5ct_immage_soundlife_combined2_sample_distribution.svg",	

rule revision:
	input:
		"plots/internal_data_prep/all_umap.svg",
		"plots/internal_data_prep/all_umap_split_by_dataset.svg",
		"plots/internal_deswan/allexp5ct_deswan_qc_plots.svg",
		expand("plots/internal_clustering/allexp5ct_deswan_deg_pseudobulk_data_limma_{gender}_pca.svg", gender=["both","female","male"]),
		expand("plots/internal_celltype_prop/{level}_propeller_{gender}_vlnplot.svg", level=["l1","l2"], gender=["both","female","male"]),
		"plots/internal_celltype_prop/celltype_prop_data_raw_scatter.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.svg",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.svg",
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3.pdf", db=["GO","Reactome"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_specific_clusters_fa_{db}_network.pdf", db=["GO","Reactome"]),
		"plots/internal_deswan/allexp5ct_deswan_removebatcheffect_q_res_all.svg",
	
		"plots/internal_deswan/allexp5ct_deswan_res_downsample.svg",
		"plots/internal_deswan/allexp5ct_deswan_q_res_by_gender2_level2.pdf",	
		expand("tables/internal_clock/allexp5ct_deswan_deg_lodo_2dmlp_{gender}_model_res.csv", gender=["both","female","male"]),
		expand("tables/external_clock/allexp5ct_deswan_deg_lodo_2dmlp_{gender}_res.csv", gender=["both","female","male"]),
		"plots/internal_clock/allexp5ct_sample_distribution_male_5cv.svg",	
		"plots/internal_clock/allexp5ct_sample_distribution_male_lodo.svg",	

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
include: "rules/external_immage_data_prep.smk"
include: "rules/external_immage_analysis.smk"
include: "rules/external_soundlife_data_prep.smk"
