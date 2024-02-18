#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_6_HCC;
done