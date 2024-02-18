#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_16_HCC;
done