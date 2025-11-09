rule test_external_deswan_enet_cv:
	input:
		"data/external_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model.joblib",
	output:
		"tables/external_clock/{mode}_deswan_deg_enet_{gender}_cv_predictions.csv",
	conda: "../env/enet.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_clock/test_deswan_enet_cv.py"

rule test_external_both_deswan_xgboost_cv:
	input:
		"data/external_clock/{mode}_deswan_deg_pseudobulk_both_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_params.json",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/external_clock/{mode}_deswan_deg_xgboost_both_cv_predictions.csv",
		"tables/external_clock/{mode}_deswan_deg_xgboost_both_cv_shap.csv",
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_clock/test_deswan_xgboost_cv.py"

rule test_external_gender_deswan_xgboost_cv:
	input:
		"data/external_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_params.json",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/external_clock/{mode}_deswan_deg_xgboost_{gender}_cv_predictions.csv",
		"tables/external_clock/{mode}_deswan_deg_xgboost_{gender}_cv_shap.csv",
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_clock/test_deswan_xgboost_cv.py"

rule test_external_both_deswan_2dmlp_cv:
	input:
		"data/external_clock/{mode}_deswan_deg_pseudobulk_both_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_scalers.joblib",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/external_clock/{mode}_deswan_deg_2dmlp_both_cv_predictions.csv",
		"tables/external_clock/{mode}_deswan_deg_2dmlp_both_cv_shap.csv",
	log:
		"log/external_clock/test_2dmlp_cv_{mode}_deswan_deg_both.log"
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_clock/test_deswan_2dmlp_cv.py"

rule test_external_gender_deswan_2dmlp_cv:
	input:
		"data/external_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_scalers.joblib",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/external_clock/{mode}_deswan_deg_2dmlp_{gender}_cv_predictions.csv",
		"tables/external_clock/{mode}_deswan_deg_2dmlp_{gender}_cv_shap.csv",
	log:
		"log/external_clock/test_2dmlp_cv_{mode}_deswan_deg_{gender}.log"
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/external_clock/test_deswan_2dmlp_cv.py"

rule compare_external_validation_deswan_all:
        input:
                pred_b1="tables/external_clock/{mode}_deswan_deg_2dmlp_both_cv_predictions.csv",
                pred_f1="tables/external_clock/{mode}_deswan_deg_2dmlp_female_cv_predictions.csv",
                pred_m1="tables/external_clock/{mode}_deswan_deg_2dmlp_male_cv_predictions.csv",
                pred_b2="tables/external_clock/{mode}_deswan_deg_enet_both_cv_predictions.csv",
                pred_f2="tables/external_clock/{mode}_deswan_deg_enet_female_cv_predictions.csv",
                pred_m2="tables/external_clock/{mode}_deswan_deg_enet_male_cv_predictions.csv",
                pred_b3="tables/external_clock/{mode}_deswan_deg_xgboost_both_cv_predictions.csv",
                pred_f3="tables/external_clock/{mode}_deswan_deg_xgboost_female_cv_predictions.csv",
                pred_m3="tables/external_clock/{mode}_deswan_deg_xgboost_male_cv_predictions.csv",
        output:
                plot1="plots/external_clock/{mode}_deswan_deg_model_comparison_metrics.pdf",
                plot2="plots/external_clock/{mode}_deswan_deg_model_comparison_scatter.pdf",
                plot3="plots/external_clock/{mode}_deswan_deg_model_comparison_disease.pdf",
		sum="tables/external_clock/{mode}_deswan_deg_model_comparison_all_summary.csv",
		fold="tables/external_clock/{mode}_deswan_deg_model_comparison_all_folds.csv",
        conda:  "../env/internal_downstream.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "normal"
        script:
                "../scripts/external_clock/compare_external_validation_all.R"
