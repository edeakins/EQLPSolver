#!/bin/bash

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_17_MLIB;
done