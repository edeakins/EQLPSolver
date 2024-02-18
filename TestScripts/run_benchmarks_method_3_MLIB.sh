#!/bin/bash

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --time_limit=3600 --report_time_file=run_method_3_MLIB;
done