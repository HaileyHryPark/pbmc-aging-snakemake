rule subset_deswan_deg_pseudobulk_data_both:
	input:
		data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/{mode}_deswan_q_deg_both.csv",
	output:
		res="data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	conda: "../env/internal_clock.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/internal_deswan/subset_deswan_deg_pseudobulk_data.R"

rule subset_deswan_deg_pseudobulk_data_female:
	input:
		data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/{mode}_deswan_q_deg_female.csv",
	output:
		res="data/internal_clock/{mode}_deswan_deg_pseudobulk_female_data_all.csv",
	conda: "../env/internal_clock.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/internal_deswan/subset_deswan_deg_pseudobulk_data.R"

rule subset_deswan_deg_pseudobulk_data_male:
	input:
		data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/{mode}_deswan_q_deg_male.csv",
	output:
		res="data/internal_clock/{mode}_deswan_deg_pseudobulk_male_data_all.csv",
	conda: "../env/internal_clock.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/internal_deswan/subset_deswan_deg_pseudobulk_data.R"

rule cv_train_and_test_gender_specific_enet_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model_params.json",
	params: gender="{gender}"
	conda: "../env/enet.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "120:00:00", queue = "super"
	script:
		"../scripts/internal_clock/cv_train_and_test_gender_specific_enet.py"

rule cv_train_and_test_both_subsample_enet_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_enet_both_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_both_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_both_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_enet_both_model_params.json",
	conda: "../env/enet.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "120:00:00", queue = "super"
	script:
		"../scripts/internal_clock/cv_train_and_test_both_subsample_enet.py"

rule cv_train_and_test_gender_specific_xgboost_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_params.json",
	params: gender="{gender}"
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "30:00:00", queue = "gpu-h200"
	script:
		"../scripts/internal_clock/cv_train_and_test_gender_specific_xgboost.py"

rule cv_train_and_test_both_subsample_xgboost_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_params.json",
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "30:00:00", queue = "gpu"
	script:
		"../scripts/internal_clock/cv_train_and_test_both_subsample_xgboost.py"

rule cv_train_and_test_gender_specific_2dmlp_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_scalers.joblib",
	params: gender="{gender}"
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "99:00:00", queue = "gpu-h200"
	script:
		"../scripts/internal_clock/cv_train_and_test_gender_specific_2dmlp.py"

rule cv_train_and_test_gender_both_subsample_2dmlp_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_res.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_scalers.joblib",
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "99:00:00", queue = "gpu-h200-int"
	script:
		"../scripts/internal_clock/cv_train_and_test_both_subsample_2dmlp.py"

rule crosstest_enet_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model.joblib",
	output:
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_crosstest_prediction.csv",
	params: gender="{gender}"
	conda: "../env/enet.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clock/crosstest_deswan_enet_cv.py"

rule crosstest_xgboost_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_params.json",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_crosstest_prediction.csv",
	params: gender="{gender}"
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu"
	script:
		"../scripts/internal_clock/crosstest_deswan_xgboost_cv.py"

rule crosstest_2dmlp_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_scalers.joblib",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_crosstest_prediction.csv",
	params: gender="{gender}"
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu"
	script:
		"../scripts/internal_clock/crosstest_deswan_2dmlp_cv.py"

rule compare_internal_validation_deswan:
	input:
		pred_b="tables/internal_clock/{mode}_deswan_deg_{model}_both_model_prediction.csv",
		pred_f="tables/internal_clock/{mode}_deswan_deg_{model}_female_model_prediction.csv",
		pred_m="tables/internal_clock/{mode}_deswan_deg_{model}_male_model_prediction.csv",
	output:
		predplot="plots/internal_clock/{mode}_deswan_deg_{model}_pred_comparison.pdf",
		resplot="plots/internal_clock/{mode}_deswan_deg_{model}_metric_comparison.pdf",
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/internal_clock/compare_internal_validation.R"

rule compare_internal_validation_deswan_all:
	input:
		pred_b1="tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_prediction.csv",
		pred_f1="tables/internal_clock/{mode}_deswan_deg_2dmlp_female_model_prediction.csv",
		predct_f1="tables/internal_clock/{mode}_deswan_deg_2dmlp_female_crosstest_prediction.csv",
		pred_m1="tables/internal_clock/{mode}_deswan_deg_2dmlp_male_model_prediction.csv",
		predct_m1="tables/internal_clock/{mode}_deswan_deg_2dmlp_male_crosstest_prediction.csv",
		pred_b2="tables/internal_clock/{mode}_deswan_deg_enet_both_model_prediction.csv",
		pred_f2="tables/internal_clock/{mode}_deswan_deg_enet_female_model_prediction.csv",
		predct_f2="tables/internal_clock/{mode}_deswan_deg_enet_female_crosstest_prediction.csv",
		pred_m2="tables/internal_clock/{mode}_deswan_deg_enet_male_model_prediction.csv",
		predct_m2="tables/internal_clock/{mode}_deswan_deg_enet_male_crosstest_prediction.csv",
		pred_b3="tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_prediction.csv",
		pred_f3="tables/internal_clock/{mode}_deswan_deg_xgboost_female_model_prediction.csv",
		predct_f3="tables/internal_clock/{mode}_deswan_deg_xgboost_female_crosstest_prediction.csv",
		pred_m3="tables/internal_clock/{mode}_deswan_deg_xgboost_male_model_prediction.csv",
		predct_m3="tables/internal_clock/{mode}_deswan_deg_xgboost_male_crosstest_prediction.csv",
	output:
		plots="plots/internal_clock/{mode}_deswan_deg_model_comparison_all.pdf",
		sum="tables/internal_clock/{mode}_deswan_deg_model_comparison_all_summary.csv",
		fold="tables/internal_clock/{mode}_deswan_deg_model_comparison_all_folds.csv",
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/internal_clock/compare_internal_validation_all.R"

rule get_shap_background_2dmlp_deswan:
	input:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_female_model_prediction.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_male_model_prediction.csv",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_shap_background.csv",
	params: both_fold=1, female_fold=4, male_fold=3
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 100, walltime = "90:00:00", queue = "gpu-h200"
	script:
		"../scripts/internal_clock/get_shap_background_2dmlp_deswan.py"
	
rule compute_shap_both_final_2dmlp_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_shap_background.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_scalers.joblib",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_shap_values.csv",
	params: gender="both" 
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 100, walltime = "90:00:00", queue = "gpu"
	script:
		"../scripts/internal_clock/compute_shap_2dmlp.py"

rule compute_shap_gender_specific_final_2dmlp_deswan:
	input:
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_shap_background.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_scalers.joblib",
	output:
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_shap_values.csv",
	params: gender="{gender}" 
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 100, walltime = "90:00:00", queue = "gpu-h200"
	script:
		"../scripts/internal_clock/compute_shap_2dmlp.py"

rule plot_shap_gender_specific_final_2dmlp_deswan:
        input:
                clust_b="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
                clust_f="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		clust_m="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
                pred_b="tables/external_clock/{mode}_deswan_deg_2dmlp_both_model_prediction_corrected_all.csv",
                pred_f="tables/external_clock/{mode}_deswan_deg_2dmlp_female_model_prediction_corrected_all.csv",
                pred_m="tables/external_clock/{mode}_deswan_deg_2dmlp_male_model_prediction_corrected_all.csv",
                shap_b="tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_shap_values.csv",
                shap_f="tables/internal_clock/{mode}_deswan_deg_2dmlp_female_model_shap_values.csv",
                shap_m="tables/internal_clock/{mode}_deswan_deg_2dmlp_male_model_shap_values.csv",
        output:
                pcaplot="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_pcaplots.svg",
                pcaplot_c="plots/internal_clock/{mode}_deswan_deg_2dmlp_corrected_shap_pcaplots.svg",
                corrplot_f="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_age_corrplots_female.svg",
                corrplot_m="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_age_corrplots_male.svg",
                vlnplot_f="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_vlnplots_female.svg",
		vlnplot_m="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_vlnplots_male.svg",
                vlnplot2_f="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_ag_vlnplots_female.svg",
                vlnplot2_m="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_ag_vlnplots_male.svg",
                vlnplot3_f="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_ag_accel_vlnplots_female.svg",
                vlnplot3_m="plots/internal_clock/{mode}_deswan_deg_2dmlp_shap_clust_ag_accel_vlnplots_male.svg",
        conda: "../env/internal_downstream.yaml"
        threads: 1
        resources: ngpus = 1, mem_gb = 20, walltime = "10:00:00", queue = "gpu"
        script:
                "../scripts/internal_clock/plot_shap_gender_specific_final_2dmlp_deswan.R"


