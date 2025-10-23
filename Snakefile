
configfile: "config.yaml"


rule internal_prep:
	input:
		"tables/internal_data_prep/internal_data_initial_metadata_table.csv",
		"tables/internal_data_prep/internal_data_final_metadata_table.csv",
		expand("data/internal_data_prep/split_h5ad_by_donor_log_{dataset}.txt", dataset=config["internal_ds"]),
		"tables/internal_data_prep/final_data_included_summary.csv",
		"data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
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

rule internal_clock:
	input:
		expand("plots/internal_clock/allexp5ct_deswan_deg_{model}_pred_comparison.pdf", model=["enet","xgboost","2dmlp"]),
		
rule internal_clustering:
	input:
		#expand("tables/internal_clustering/allexp5ct_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv", gender=["both","female","male"]),
		expand("plots/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_fa_{gender}_res.csv", gender=["both","female","male"]),
		expand("tables/internal_clustering/allexp5ct_deswan_deg_mfuzz_merged_clusters_disgenet2r_{gender}_res.csv", gender=["both","female","male"]),

##### load rules #####

include: "rules/r_package_install.smk"
include: "rules/internal_data_prep.smk"
include: "rules/internal_pseudobulk.smk"
include: "rules/internal_deswan.smk"
include: "rules/internal_clustering.smk"
include: "rules/internal_clock.smk"
