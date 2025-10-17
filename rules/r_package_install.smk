rule install_github_r_package:
	input:
		"resources/r_package"
	output:
		directory("resources/r_package/{pkg}")
	params:
		path=lambda w: config["r_package"][w.pkg]["path"]
	log:
		"log/preprocessing/install_{pkg}.log"
	localrule: True
	conda: "../env/r_package_install.yaml"
	script:
		"../scripts/install_r_package.R"

	


