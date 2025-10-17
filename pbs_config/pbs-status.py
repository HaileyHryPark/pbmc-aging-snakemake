#!/usr/bin/env python3
# Adapted from https://raw.githubusercontent.com/Snakemake-Profiles/pbs-torque/master/%7B%7Bcookiecutter.profile_name%7D%7D/pbs-status.py

import re, sys
import subprocess

jobid = sys.argv[1]

try:
    res = subprocess.run(
        "qstat -f -x {}".format(jobid), check=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, shell=True)
    job_state = re.search(
        "(?<=job_state = )[A-Z]", res.stdout.decode()).group(0)

    if job_state == "F":  # F means finished
        # Had to add r in front of the string for a `SyntaxWarning: invalid escape sequence`
        # https://stackoverflow.com/questions/57645829/why-am-i-getting-a-syntaxwarning-invalid-escape-sequence-s-warning
        exit_status = re.search(r"(?<=Exit_status = )\d+", res.stdout.decode())
        print(
            # Finished, but exit_status not yet assigned; associated with jobs
            # held for a long time, and qsub gives up trying to run it.
            "running" if exit_status is None else
            "success" if exit_status.group(0) == "0" else "failed")
    else:
        print("running")

except (subprocess.CalledProcessError, IndexError, KeyboardInterrupt) as e:
    print("failed")
