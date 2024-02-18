#!/bin/bash

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=simplex --time_limit=3600 --report_time_file=run_method_1_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=simplex --time_limit=3600 --report_time_file=run_method_1_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm --time_limit=3600 --report_time_file=run_method_2_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm --time_limit=3600 --report_time_file=run_method_2_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --time_limit=3600 --report_time_file=run_method_3_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual --time_limit=3600 --report_time_file=run_method_3_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm --time_limit=3600 --report_time_file=run_method_4_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm --time_limit=3600 --report_time_file=run_method_4_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_5_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_5_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_6_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_6_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_7_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_7_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=impcross_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_8_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=impcross_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_8_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_9_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_standard_crossover --time_limit=3600 --report_time_file=run_method_9_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_10_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_primal_crossover --time_limit=3600 --report_time_file=run_method_10_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual_no_iter --time_limit=3600 --report_time_file=run_method_11_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_dual_no_iter --time_limit=3600 --report_time_file=run_method_11_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm_no_iter --time_limit=3600 --report_time_file=run_method_12_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=orbital_crossover_ipm_no_iter --time_limit=3600 --report_time_file=run_method_12_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_13_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_13_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_primal_crossover_iter --time_limit=3600 --report_time_file=run_method_14_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=dual_alp_primal_crossover_iter --time_limit=3600 --report_time_file=run_method_14_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_15_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_15_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_primal_crossover_iter --time_limit=3600 --report_time_file=run_method_16_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipmcross_alp_primal_crossover_iter --time_limit=3600 --report_time_file=run_method_16_MLIB;
done

# HS-COV-COD Instances
for f in ../HS-COV-COD/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_17_HCC;
done

# MIPLIB Instances
for f in ../MIPLIB-2017/*; do
   ../HiGHS-1-2-1/build/bin/highs "$f" --solver=ipm_alp_standard_crossover_iter --time_limit=3600 --report_time_file=run_method_17_MLIB;
done

