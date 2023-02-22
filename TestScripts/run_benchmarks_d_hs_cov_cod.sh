#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --report_time_file=HiGHS_DUAL_HS-COV-COD_TIMES --time_limit=3600;
done