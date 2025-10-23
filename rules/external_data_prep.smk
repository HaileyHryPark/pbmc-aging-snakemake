rule prep_marina_data:
        input:
                "data/internal_data_prep/all_pbmcs_rna.h5ad",
                "data/internal_data_prep/all_pbmcs_metadata.csv",
        output:
                "data/internal_data_prep/marina.h5ad",
	conda: "../env/internal_data_prep_py.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
        script:
                "../scripts/internal_data_prep/prep_marina_data.py"

rule extract_metadata:
        input:
                "data/internal_data_prep/{dataset}.h5ad",
        output:
                "data/internal_data_prep/{dataset}_qced.h5ad",
                "data/internal_data_prep/{dataset}_metadata.csv",
                "data/internal_data_prep/{dataset}_ensembl_to_symbol.csv",
	params: dataset="{dataset}"
	conda: "../env/internal_data_prep_py.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 200, walltime = "20:00:00", queue = "normal"
        script:
                "../scripts/internal_data_prep/extract_metadata.py"

rule write_initial_metadata_table:
        input:
                onek1k="data/internal_data_prep/onek1k_metadata.csv",
                aida="data/internal_data_prep/aida_metadata.csv",
                marina="data/internal_data_prep/marina_metadata.csv",
                perez="data/internal_data_prep/perez_metadata.csv",
        output:
                table="tables/internal_data_prep/internal_data_initial_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
        script:
                "../scripts/internal_data_prep/write_metadata_table.R"

rule filter_data:
        input:
                "data/internal_data_prep/{dataset}_qced.h5ad",
                "tables/internal_data_prep/internal_data_initial_metadata_table.csv",
        output:
                "data/internal_data_prep/{dataset}_filtered.h5ad",
                "data/internal_data_prep/{dataset}_filtered_metadata.csv",
	conda: "../env/internal_data_prep_py.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
        script:
                "../scripts/internal_data_prep/filter_data.py"

rule write_final_metadata_table:
        input:
                onek1k="data/internal_data_prep/onek1k_filtered_metadata.csv",
                aida="data/internal_data_prep/aida_filtered_metadata.csv",
                marina="data/internal_data_prep/marina_filtered_metadata.csv",
                perez="data/internal_data_prep/perez_filtered_metadata.csv",
        output:
                table="tables/internal_data_prep/internal_data_final_metadata_table.csv",
	conda: "../env/internal_data_prep.yaml"
        threads: 1
        resources: ngpus = 0, mem_gb = 200, walltime = "02:00:00", queue = "normal"
        script:
                "../scripts/internal_data_prep/write_metadata_table.R"

#rule split_h5ad_by_donor:
#
## manual_conversion_h5ad_to_seu
#
#rule run_azimuth:

