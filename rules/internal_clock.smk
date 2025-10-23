rule subset_deswan_deg_pseudobulk_data_both:
        input:
                data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
                deg="tables/internal_deswan/{mode}_deswan_q_deg_both.csv",
        output:
                res="data/internal_clock/{mode}_deswan_deg_pseudobulk_data_all.csv",
        conda: "../env/internal_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "normal"
        script:
                "../scripts/internal_deswan/subset_deswan_deg_pseudobulk_data.R"

rule subset_deswan_deg_pseudobulk_data_female:
        input:
                data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
                deg="tables/internal_deswan/{mode}_deswan_q_deg_female.csv",
        output:
                res="data/internal_clock/{mode}_deswan_deg_pseudobulk_female_data_all.csv",
        conda: "../env/internal_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "normal"
        script:
                "../scripts/internal_deswan/subset_deswan_deg_pseudobulk_data.R"

rule subset_deswan_deg_pseudobulk_data_male:
        input:
                data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
                deg="tables/internal_deswan/{mode}_deswan_q_deg_male.csv",
        output:
                res="data/internal_clock/{mode}_deswan_deg_pseudobulk_male_data_all.csv",
        conda: "../env/internal_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 150, walltime = "02:00:00", queue = "normal"
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
	resources: ngpus = 0, mem_gb = 120, walltime = "120:00:00", queue = "normal"
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
	resources: ngpus = 0, mem_gb = 120, walltime = "120:00:00", queue = "normal"
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
	resources: ngpus = 1, mem_gb = 120, walltime = "30:00:00", queue = "normal"
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
	resources: ngpus = 1, mem_gb = 120, walltime = "30:00:00", queue = "normal"
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
	resources: ngpus = 1, mem_gb = 120, walltime = "99:00:00", queue = "normal"
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
	resources: ngpus = 1, mem_gb = 120, walltime = "99:00:00", queue = "normal"
	script:
		"../scripts/internal_clock/cv_train_and_test_both_subsample_2dmlp.py"

rule compare_internal_validation_deswan:
	input:
		pred_b="tables/internal_clock/{mode}_deswan_deg_{model}_both_model_prediction.csv",
		pred_f="tables/internal_clock/{mode}_deswan_deg_{model}_female_model_prediction.csv",
		pred_m="tables/internal_clock/{mode}_deswan_deg_{model}_male_model_prediction.csv",
	output:
		predplot="plots/internal_clock/{mode}_deswan_deg_{model}_pred_comparison.pdf",
	conda:	"../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "05:00:00", queue = "normal"
	script:
		"../scripts/internal_clock/compare_internal_validation.R"

