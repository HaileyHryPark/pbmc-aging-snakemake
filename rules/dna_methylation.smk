rule qc_before_filter: 
	input:
		meta="data/dna_methylation/samples.csv",
	output:
		metadata="data/dna_methylation/nsphs_metadata.csv",
		qcplot="plots/dna_methylation/nsphs_raw_qc.pdf",
		pcaplot="plots/dna_methylation/nsphs_norm_pca.pdf",
		rgset="data/dna_methylation/nsphs_rgset.rds",
		mset="data/dna_methylation/nsphs_norm.rds",
	conda: "../env/dna_methylation_qc.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/qc_before_filter.R"

rule plot_sample_distribution_dnam_nsphs:
        input:
                metadata="data/dna_methylation/nsphs_metadata.csv",
        output:
                plot="plots/dna_methylation/nsphs_sample_distribution_age_gender.pdf",
        conda: "../env/dna_methylation_qc.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
        script:
                "../scripts/dna_methylation/plot_sample_distribution.R"

rule filter_probes: 
	input:
		rgset="data/dna_methylation/nsphs_rgset.rds",
		mset="data/dna_methylation/nsphs_norm.rds",
		nsref="resources/dna_methylation/48639-non-specific-probes-Illumina450k.csv",
		mmref="resources/dna_methylation/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt",
	output:
		filtered="data/dna_methylation/nsphs_filtered.rds",
		probes="data/dna_methylation/nsphs_probes.csv",
	conda: "../env/dna_methylation_qc.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/filter_probes.R"

rule get_beta_mvalue:
	input:
		mset="data/dna_methylation/nsphs_filtered.rds",
		probes="data/dna_methylation/nsphs_probes.csv",
	output:
		beta="tables/dna_methylation/nsphs_beta.csv",
		mvalue="tables/dna_methylation/nsphs_mvalue.csv",
	conda: "../env/dna_methylation_beta.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/get_beta_mvalue.R"

rule subset_data_by_gender:
	input:
		beta="tables/dna_methylation/nsphs_beta.csv",
		mvalue="tables/dna_methylation/nsphs_mvalue.csv",
	output:
		male_b="tables/dna_methylation/nsphs_beta.csv",
		mvalue="tables/dna_methylation/nsphs_mvalue.csv",
	conda: "../env/dna_methylation_beta.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/get_beta_mvalue.R"

rule run_limma_dmp:
	input:
		metadata="data/dna_methylation/nsphs_metadata.csv",
		mvalue="tables/dna_methylation/nsphs_mvalue.csv",
	output:
		res1="tables/dna_methylation/nsphs_limma_dmp_age_res.csv",
		res2="tables/dna_methylation/nsphs_limma_dmp_agegroup_res.csv",
	conda: "../env/cluster_score.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/run_limma_dmp.R"

rule annotate_limma_dmp_res: 
	input:
		limma="tables/dna_methylation/nsphs_limma_dmp_{design}_res.csv",
	output:
		annotated="tables/dna_methylation/nsphs_limma_dmp_{design}_res_annotated.csv",
	conda: "../env/dna_methylation_qc.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/annotate_limma_dmp_res.R"

rule run_limma_dmrcate_gender_aware:
	input:
		metadata="data/dna_methylation/nsphs_metadata.csv",
		mvalue="tables/dna_methylation/nsphs_mvalue.csv",
	output:
		res_both="tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_both.rds",
		res_female="tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_female.rds",
		res_male="tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_male.rds",
	conda: "../env/dna_methylation_dm.yaml"
	threads: 1
	resources: ngpus = 0, mem_gb = 120, walltime = "20:00:00", queue = "normal"
	script:
		"../scripts/dna_methylation/run_limma_dmrcate_gender_aware.R"

rule plot_limma_dmrcate_sig_gender_interaction:
        input:
                res="tables/dna_methylation/nsphs_limma_dmrcate_agegroup_res_both.rds",
                metadata="data/dna_methylation/nsphs_metadata.csv",
                beta="tables/dna_methylation/nsphs_beta.csv",
        output:
                pie="plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_chromosome_pie.pdf",
                plot1="plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_plots.pdf",
                plot2="plots/dna_methylation/nsphs_limma_dmrcate_sig_gender_interaction_line_vln_plots.pdf",
        conda: "../env/dna_methylation_dm.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 50, walltime = "08:00:00", queue = "normal"
        script:
                "../scripts/dna_methylation/plot_limma_dmrcate_sig_gender_interaction.R"
	
