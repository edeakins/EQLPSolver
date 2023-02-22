#!/bin/bash

# MIPLIB 2017 Instnaces
for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm --report_time_file=OC_IPM_MIPLIB-2017_TIMES --time_limit=3600;
done

for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --report_time_file=OC_DUAL_MIPLIB-2017_TIMES --time_limit=3600;
done

for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_aggregate --report_time_file=HiGHS_IPM_ALP_MIPLIB-2017_TIMES --time_limit=3600;
done

for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm --report_time_file=HiGHS_IPM_MIPLIB-2017_TIMES --time_limit=3600;
done

for f in ../MIPLIB-2017/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --report_time_file=HiGHS_DUAL_MIPLIB-2017_TIMES --time_limit=3600;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm --report_time_file=OC_IPM_HS-COV-COD_TIMES --time_limit=3600;
done

for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --report_time_file=OC_DUAL_HS-COV-COD_TIMES --time_limit=3600;
done

for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_aggregate --report_time_file=HiGHS_IPM_ALP_HS-COV-COD_TIMES --time_limit=3600;
done

for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm --report_time_file=HiGHS_IPM_HS-COV-COD_TIMES --time_limit=3600;
done

for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --report_time_file=HiGHS_DUAL_HS-COV-COD_TIMES --time_limit=3600;
done

