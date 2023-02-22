#!/bin/bash

# MIPLIB 2017 Instnaces
for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --report_time_file=OC_DUAL_MIPLIB-2017_TIMES --time_limit=3600;
done