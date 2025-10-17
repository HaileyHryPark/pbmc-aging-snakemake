
configfile: "config.yaml"


rule all:
	input:


##### load rules #####

include: "rules/r_package_install.smk"
include: "rules/rule.smk"
