#!/bin/bash
for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/Release/bin/highs "$f" --solver=orbital_crossover_dual --report_time_file=OC_DUAL_HS-COV-COD_TIMES_ONE_LIFT --time_limit=3600;
done