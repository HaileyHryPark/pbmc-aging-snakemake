# snakemake-project-template

Initial conda environment for snakemake configuration at env/initial_env.yaml

Profile set up for hpc usage in pbs_config/
 
Can run these commands to run snakemake in hpc:
```
casmk
```
```
smk
```

Aliases:
```
alias casmk='conda activate smk'
alias smk='snakemake --profile ./pbs_config --use-conda --latency-wait 60'
```

