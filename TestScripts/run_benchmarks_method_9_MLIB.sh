#!/bin/bash

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=impcross_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_9_MLIB;
done