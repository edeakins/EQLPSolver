#!/bin/bash

for f in ../mittleman_symmetric_1/*; do
	../DHiGHS/build/bin/highs "$f" --aggregate=on --time_file=MIPLIB2010_OC_times.csv --reduction_file=MIPLIB2010_reductions.csv --time_limit=3600;
	../DHiGHS/build/bin/highs "$f" --time_file=MIPLIB2010_HiGHS_times.csv --time_limit=3600;
done

for f in ../symmetric_1/*; do
	../DHiGHS/build/bin/highs "$f" --aggregate=on --time_file=HS-COV-COD_OC_times.csv --reduction_file=HS-COV-COD_reductions.csv --time_limit=3600;
	../DHiGHS/build/bin/highs "$f" --time_file=HS-COV-COD_HiGHS_times.csv --time_limit=3600;
done

