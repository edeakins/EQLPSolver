#!/bin/bash
for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm --report_time_file=HiGHS_IPM_HS-COV-COD_TIMES --time_limit=3600;
done