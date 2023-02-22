#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm --report_time_file=OC_IPM_HS-COV-COD_TIMES --time_limit=3600;
done