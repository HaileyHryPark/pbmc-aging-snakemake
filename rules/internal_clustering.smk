rule run_limma_deswan_pseudobulk_data:
	input:
		data="data/internal_pseudobulk/{mode}_pseudobulk_data_all.csv",
		deg="tables/internal_deswan/{mode}_deswan_q_deg_{gender}.csv",
	output:
		res="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}.csv",
	params: gender="{gender}"
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 80, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_limma_removebatcheffect_pseudobulk_data.R"

rule fit_gender_loess_zscore_limma_deswan_deg_subset:
	input:
		limma="data/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}.csv",
	output:
		zscaled="tables/internal_clustering/{mode}_deswan_deg_pseudobulk_data_limma_{gender}_zscaled.csv",
		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_{gender}.csv",
		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_{gender}.csv",
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 100, walltime = "80:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/fit_gender_loess_zscore_limma_deswan_deg_subset.R"

rule estimate_params_mfuzz_loess_fitted:
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.pdf",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/estimate_params_mfuzz_loess_fitted.R"

rule cluster_mfuzz_loess_fitted_multirun: 
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
	output:
		cluster_counts="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_cluster_num_{gender}.csv",
		all_runs_rds="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_res_{gender}.rds",
	conda: "../env/internal_clustering2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/cluster_mfuzz_loess_fitted_multirun.R"

rule cluster_mfuzz_loess_fitted_final_run: 
	input:
		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
	output:
		final_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_assignment_{gender}.csv",
		final_centers="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_centers_{gender}.csv",
	conda: "../env/internal_clustering2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/cluster_mfuzz_loess_fitted_final_run.R"

rule plot_mfuzz_loess_fitted_cluster_runs: 
	input:
		cluster_counts="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_cluster_num_{gender}.csv",
		all_runs_rds="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_res_{gender}.rds",
		final_centers="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_centers_{gender}.csv",
	output:
		plot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_cluster_num_{gender}.svg",
		plot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_trajectories_{gender}.svg",
		plot3="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_multirun_cluster_trajectories_{gender}.svg",
	conda: "../env/internal_clustering2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_loess_fitted_cluster_runs.R"

rule annotate_mfuzz_loess_fitted_final_clusters: 
	input:
		final_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_finalrun_cluster_assignment_{gender}.csv",
	output:
		annotated="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	params: gender="{gender}"
	conda: "../env/internal_clustering2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/annotate_mfuzz_loess_fitted_final_clusters.R"

#rule cluster_mfuzz_loess_fitted_final:
#	input:
#		mfuzz_mat="tables/internal_clustering/{mode}_deswan_deg_loess_mfuzz_mat_{gender}.txt",
#		table="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_num_{gender}.csv",
#	output:
#		cnumplot="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_elbowplot_{gender}.pdf",
#		corrplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_initial_{gender}.pdf",
#		mfuzzplot1="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_initial_{gender}.pdf",
#		corrplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_corrplot_merged_{gender}.pdf",
#		mfuzzplot2="plots/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_plot_merged_{gender}.pdf",
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv",
#	conda: "../env/internal_clustering.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
#	script:
#		"../scripts/internal_clustering/cluster_mfuzz_loess_fitted_final.R"
#
#rule plot_mfuzz_merged_clusters_loess:
#	input:
#		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_{gender}.csv",
#		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_{gender}.csv",
#		var_cluster_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}.csv",
#	output:
#		annotated="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
#		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_loess_{gender}.pdf",
#	params: gender="{gender}"
#	conda: "../env/internal_clustering.yaml"
#	threads: 1
#	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
#	script:
#		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_loess.R"

rule plot_mfuzz_merged_clusters_loess_final:
	input:
		fit_res="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_{gender}.csv",
		span_res="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_{gender}.csv",
		annotated="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_{gender}_annotated.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_loess_{gender}_final.svg",
	params: gender="{gender}"
	conda: "../env/final_plots.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_loess_final.R"

rule plot_mfuzz_merged_clusters_info:
	input:
		both="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		flowplot="plots/internal_clustering/{mode}_deswan_deg_gender_flowplot_by_merged_clusters.pdf",
		venn="plots/internal_clustering/{mode}_deswan_deg_gender_venn_by_merged_clusters.pdf",
		pie="plots/internal_clustering/{mode}_deswan_cluster_pie_by_celltype.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_info.R"

rule run_mfuzz_merged_clusters_fa_all:
	input:
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res_no_ribo.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_fa_all.R"

rule plot_mfuzz_merged_clusters_fa_all:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv",
		annot="data/internal_clustering/top_fa_terms.csv",
	output:
		go="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_go_res_plots.pdf",
		annotgo="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_go_top_res_annotated.csv",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa_all.R"

rule plot_mfuzz_merged_clusters_fa_score_all:
	input:
		data="data/internal_pseudobulk/full5ct_pseudobulk_data_all.csv",
		annotgo="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_go_top_res_annotated.csv",
	output:
		scores="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_top_fa_scores.csv",
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_top_fa_score_plots.pdf",
	conda: "../env/cluster_score.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "99:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa_singscore_all.R"

rule plot_mfuzz_merged_clusters_fa_network:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv",
	output:
		all="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res_network_plots.pdf",
		go="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_go_res_network_plots.pdf",
		wp="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_wp_res_network_plots.pdf",
		r="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_r_res_network_plots.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa_network.R"

rule plot_mfuzz_merged_clusters_fa_network2:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv",
		annot="data/internal_clustering/top_fa_terms.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res_network_plots2.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa_network2.R"

rule plot_mfuzz_merged_clusters_fa_network3:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_all_res.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3.pdf",
		plot_fli="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3_fli.pdf",
                plot_mef="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_{db}_res_network3_mef.pdf",
	params: db="{db}"
	conda: "../env/internal_downstream2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa_network3.R"

rule run_mfuzz_merged_clusters_msigdb_fa_all:
	input:
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_all_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_all_res_no_ribo.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_msigdb_fa_all.R"

rule plot_mfuzz_merged_clusters_msigdb_fa_all:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_all_res.csv",
	output:
		chr="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_msigdb_fa_chr_res_plots.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_msigdb_fa_all.R"

rule run_mfuzz_merged_clusters_disgenet2r_all:
	input:
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_all_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_disgenet2r_all_res_no_ribo.csv",
	singularity: "/apps/singularity/rstudio-4.5.0_ExtPack_NOV102025.sif"
	#conda: "../env/azimuth.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_disgenet2r_all.R"

rule plot_mfuzz_merged_clusters_fa:
	input:
		both1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_both_res.csv",
		both2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_both_res_no_ribo.csv",
		female1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_female_res.csv",
		female2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_female_res_no_ribo.csv",
		male1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_male_res.csv",
		male2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_male_res_no_ribo.csv",
	output:
		plot1="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_plot.pdf",
		plot2="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_fa_plot_no_ribo.pdf",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_fa.R"

rule run_mfuzz_merged_clusters_mitocarta_fa_all:
	input:
		mitocarta="resources/MitoCarta3.0_data.rds",
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_mitocarta_fa_all_res.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_merged_clusters_mitocarta_fa_all.R"

rule plot_mfuzz_merged_clusters_mitocarta_fa_all:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_mitocarta_fa_all_res.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_merged_clusters_mitocarta_fa_res_plots.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_merged_clusters_mitocarta_fa_all.R"

rule plot_mfuzz_specific_clusters_loess:
	input:
		fit_res_f="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_female.csv",
		span_res_f="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_female.csv",
		clust_df_f="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		fit_res_m="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_male.csv",
		span_res_m="tables/internal_clustering/{mode}_deswan_deg_loess_span_res_male.csv",
		clust_df_m="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		plot1="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_loess.pdf",
		plot2="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_celltype_pie.pdf",
	conda: "../env/internal_clustering.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_specific_clusters_loess.R"

rule run_mfuzz_specific_clusters_fa_all:
	input:
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		res1="tables/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_all_res.csv",
		res2="tables/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_all_res_no_ribo.csv",
	conda: "../env/functional_annotation.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_specific_clusters_fa_all.R"

rule plot_mfuzz_specific_clusters_fa_all:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_all_res.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_all_res.pdf",
	conda: "../env/internal_downstream.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_specific_clusters_fa_all.R"

rule plot_mfuzz_specific_clusters_fa_network:
	input:
		table="tables/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_all_res.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_{db}_network.pdf",
		plot_fcimei="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_fa_{db}_network_fcimei.pdf",
	params: db="{db}"
	conda: "../env/internal_downstream2.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/plot_mfuzz_specific_clusters_fa_network.R"

rule run_mfuzz_specific_clusters_disgenet2r_all:
	input:
		both_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_both_annotated.csv",
		female_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_female_annotated.csv",
		male_df="tables/internal_clustering/{mode}_deswan_deg_loess_fitted_mfuzz_cluster_assignment_male_annotated.csv",
	output:
		plot="plots/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_disgenet2r_all_res.pdf",
		res="tables/internal_clustering/{mode}_deswan_deg_mfuzz_specific_clusters_disgenet2r_all_res.csv",
	conda: "../env/azimuth.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 50, walltime = "20:00:00", queue = "super"
	script:
		"../scripts/internal_clustering/run_mfuzz_specific_clusters_disgenet2r_all.R"

