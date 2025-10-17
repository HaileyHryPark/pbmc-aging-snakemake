rule test: 
	input:
		data=""
	output:
		obj=""
		plot=report(
			"plots/",
			caption="report/report.rst",
			category="",
		),
		table=report(
			"tables/",
			caption="report/report.rst"
			category="",
		)
	params:
		min_res=lambda w: config[w.dataset]["clustering"]["min_res"],
	log:
		"log/"
	conda: "../env/file.yaml"
	threads: 20
	resources: ngpus = 1, mem_gb = 100, walltime = "24:00:00", queue = "gpu"
	script: 
		"../scripts/test.R"




