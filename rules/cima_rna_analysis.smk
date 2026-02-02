rule run_deswan_cima_rna:
	input:
		data="data/cima_rna_data_prep/{mode}_pseudobulk_data_cima_rna.csv",
	output:
		res="tables/cima_rna_analysis/{mode}_deswan_res_cima_rna.rds",
		coef="tables/cima_rna_analysis/{mode}_deswan_coef_res_cima_rna.csv",
		p="tables/cima_rna_analysis/{mode}_deswan_p_res_cima_rna.csv",
		q="tables/cima_rna_analysis/{mode}_deswan_q_res_cima_rna.csv",
	#conda: "../env/internal_data_prep.yaml"
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "99:00:00", queue = "long"
	script:
		"../scripts/cima_rna_analysis/run_deswan.R"

rule plot_deswan_cima_rna:
	input:
		q="tables/cima_rna_analysis/{mode}_deswan_q_res_cima_rna.csv",
	output:
		plots="plots/cima_rna_analysis/{mode}_deswan_q_res_by_gender.pdf",
		plots2="plots/cima_rna_analysis/{mode}_deswan_q_res_by_gender2.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/plot_deswan.R"

rule get_deswan_deg_cima_rna:
	input:
		q="tables/cima_rna_analysis/{mode}_deswan_q_res_cima_rna.csv",
	output:
		degb="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_both.csv",
		degf="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_female.csv",
		degm="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_male.csv",
	conda: "../env/internal_deswan.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/get_deswan_deg.R"	 

rule plot_deswan_deg_venn_cima_rna:
	input:
		degb="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_both.csv",
		degf="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_female.csv",
		degm="tables/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_male.csv",
	output:
		venn="plots/cima_rna_analysis/{mode}_deswan_q_deg_cima_rna_venn.pdf"
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "02:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/plot_deswan_deg_venn.R"	 

rule test_deswan_enet_cv_cima_rna:
	input:
		"data/cima_rna_data_prep/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_enet_{gender}_model.joblib",
	output:
		"tables/cima_rna_analysis/{mode}_deswan_deg_enet_{gender}_cv_predictions.csv",
	conda: "../env/enet.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/test_deswan_enet_cv.py"

rule test_both_deswan_xgboost_cv_cima_rna:
	input:
		"data/cima_rna_data_prep/{mode}_deswan_deg_pseudobulk_both_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_both_model_params.json",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
	output:
		"tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_both_cv_predictions.csv",
		"tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_both_cv_shap.csv",
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu-h200-int"
	script:
		"../scripts/cima_rna_analysis/test_deswan_xgboost_cv.py"

rule test_gender_deswan_xgboost_cv_cima_rna:
	input:
		"data/cima_rna_data_prep/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model.joblib",
		"tables/internal_clock/{mode}_deswan_deg_xgboost_{gender}_model_params.json",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
	output:
		"tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_{gender}_cv_predictions.csv",
		"tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_{gender}_cv_shap.csv",
	conda: "../env/xgboost.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu-h200-int"
	script:
		"../scripts/cima_rna_analysis/test_deswan_xgboost_cv.py"

rule test_both_deswan_2dmlp_cv_cima_rna:
	input:
		"data/cima_rna_data_prep/{mode}_deswan_deg_pseudobulk_both_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_both_model_scalers.joblib",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_shap_background.csv",
	output:
		"tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_both_cv_predictions.csv",
		"tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_both_cv_shap.csv",
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu"
	script:
		"../scripts/cima_rna_analysis/test_deswan_2dmlp_cv.py"

rule test_gender_deswan_2dmlp_cv_cima_rna:
	input:
		"data/cima_rna_data_prep/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model.pt",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_params.json",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_scalers.joblib",
		"data/internal_clock/{mode}_deswan_deg_pseudobulk_{gender}_data_all.csv",
		"tables/internal_clock/{mode}_deswan_deg_2dmlp_shap_background.csv",
	output:
		"tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_{gender}_cv_predictions.csv",
		"tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_{gender}_cv_shap.csv",
	conda: "../env/mlp.yaml"
	threads: 1
	resources: ngpus = 1, mem_gb = 120, walltime = "20:00:00", queue = "gpu-h200"
	script:
		"../scripts/cima_rna_analysis/test_deswan_2dmlp_cv.py"

rule compare_validation_deswan_all_cima_rna:
	input:
		pred_b1="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_both_cv_predictions.csv",
		pred_f1="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_female_cv_predictions.csv",
		pred_m1="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_male_cv_predictions.csv",
		pred_b2="tables/cima_rna_analysis/{mode}_deswan_deg_enet_both_cv_predictions.csv",
		pred_f2="tables/cima_rna_analysis/{mode}_deswan_deg_enet_female_cv_predictions.csv",
		pred_m2="tables/cima_rna_analysis/{mode}_deswan_deg_enet_male_cv_predictions.csv",
		pred_b3="tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_both_cv_predictions.csv",
		pred_f3="tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_female_cv_predictions.csv",
		pred_m3="tables/cima_rna_analysis/{mode}_deswan_deg_xgboost_male_cv_predictions.csv",
	output:
		plot1="plots/cima_rna_analysis/{mode}_deswan_deg_model_comparison_metrics.pdf",
		plot2="plots/cima_rna_analysis/{mode}_deswan_deg_model_comparison_scatter.pdf",
		sum="tables/cima_rna_analysis/{mode}_deswan_deg_model_comparison_all_summary.csv",
		fold="tables/cima_rna_analysis/{mode}_deswan_deg_model_comparison_all_folds.csv",
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/compare_validation_all.R"

rule correct_2dmlp_predicted_age_cima_rna:
	input:
		pred_i="tables/internal_clock/{mode}_deswan_deg_2dmlp_{gender}_model_prediction.csv",
		pred_e="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_{gender}_cv_predictions.csv",
	output:
		corrected="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_{gender}_model_prediction_corrected_all.csv",
	params: gender="{gender}"
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/correct_2dmlp_predicted_age.R"

rule plot_corrected_2dmlp_predicted_age_cima_rna:
	input:
		both="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_both_model_prediction_corrected_all.csv",
		female="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_female_model_prediction_corrected_all.csv",
		male="tables/cima_rna_analysis/{mode}_deswan_deg_2dmlp_male_model_prediction_corrected_all.csv",
	output:
		plot1="plots/cima_rna_analysis/{mode}_deswan_deg_2dmlp_model_prediction_corrected_scatter.pdf",
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "super"
	script:
		"../scripts/cima_rna_analysis/plot_corrected_2dmlp_predicted_age.R"

