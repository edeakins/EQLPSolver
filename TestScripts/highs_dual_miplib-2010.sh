#!/bin/bash
for f in ../MIPLIB-2010/*; do
	../HiGHS-1-2-1/build/bin/highs "$f" --report_time_file=HiGHS_DUAL_MIPLIB-2010_TIMES --time_limit=3600;
done