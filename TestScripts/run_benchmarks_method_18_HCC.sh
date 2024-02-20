#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_primal_crossover_iter --time_limit=3600 --report_time_file=run_method_18_HCC;
done