
configfile: "config.yaml"


rule internal_prep:
	input:
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
		"tables/internal_data_prep/internal_data_final_metadata_table.csv",
		expand("data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt", dataset=config["internal_ds"]),
		"tables/internal_data_prep/final_data_included_summary.csv",
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
		
rule internal_clustering:
	input:
		#expand("tables/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_{gender}_res.csv", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_msigdb_fa_{gender}_res.csv", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res.csv", gender=["both","female","male"]),
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_plot.pdf",
		"plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_plot_no_ribo.pdf",

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

##### load rules #####

include: "rules/r_package_install.smk"
include: "rules/internal_data_prep.smk"
include: "rules/internal_pseudobulk.smk"
include: "rules/internal_deswan.smk"
include: "rules/internal_clustering.smk"
include: "rules/internal_clock.smk"
include: "rules/external_dis_data_prep.smk"
include: "rules/external_sc_data_prep.smk"
include: "rules/external_pseudobulk.smk"
include: "rules/external_clock.smk"
