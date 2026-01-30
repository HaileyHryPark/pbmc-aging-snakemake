rule get_cluster_score:
	input:
		data="data/{cohort}_pseudobulk/full5ct_pseudobulk_data_all.csv",
		clust_b="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		clust_f="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		clust_m="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		scores="tables/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores.csv",
	params: cohort="{cohort}"
	conda: "../env/cluster_score.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/cluster_score/get_cluster_score.R"

rule compare_age_accel_cluster_score:
	input:
		data="tables/external_clock/{mode}_deswan_deg_2dmlp_{gender}_model_prediction_corrected_all.csv",
		score="tables/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores.csv",
	output:
		boxplot="plots/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores_{gender}_age_accel_boxplot.pdf",
		radarplot="plots/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores_{gender}_age_accel_radarplot.pdf",
	params: cohort="{cohort}", gender="{gender}",
	conda: "../env/cluster_score.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/cluster_score/compare_age_accel_cluster_score.R"

rule compare_disease_cluster_score:
	input:
		score="tables/cluster_score/{mode}_deswan_deg_external_cluster_scores.csv",
	output:
		boxplot="plots/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores_{gender}_age_accel_boxplot.pdf",
		radarplot="plots/cluster_score/{mode}_deswan_deg_{cohort}_cluster_scores_{gender}_age_accel_radarplot.pdf",
	params: cohort="{cohort}", gender="{gender}",
	conda: "../env/cluster_score.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/cluster_score/compare_disease_cluster_score.R"

