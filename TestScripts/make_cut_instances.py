"""
For each solved symmetric instance, add an optimality cut:
  sum(x) >= floor(sum_x_opt) + 1
where sum_x_opt is the sum of LP solution values at optimality.
This forces the LP away from the degenerate all-at-bounds optimal,
creating non-degenerate interior solutions for testing.
Output: ../HS-COV-COD-CUT/{instance}_c.mps
"""
import highspy
import math
import gzip
import os
import sys
import tempfile
import shutil
import numpy as np

INSTANCES = [
    "cod103","codbt161","codbt162","codbt251","codbt262","codbt342",
    "codbt531","codbt532","codbt812","cov1385","cov1386","cov1387",
    "cov14107","cov14108","cov1496","cov15107","cov151110","cov15118",
    "cov161110","cov161211","cov171311","sts729",
]

os.makedirs("../HS-COV-COD-CUT", exist_ok=True)

for name in INSTANCES:
    src = f"../HS-COV-COD/{name}.mps.gz"
    dst = f"../HS-COV-COD-CUT/{name}_c.mps"

    # Decompress to temp file then solve
    tmp = tempfile.NamedTemporaryFile(suffix=".mps", delete=False)
    with gzip.open(src, 'rb') as f_in:
        shutil.copyfileobj(f_in, tmp)
    tmp.close()

    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.readModel(tmp.name)
    os.unlink(tmp.name)
    h.run()

    status = h.getModelStatus()
    if str(status) != "HighsModelStatus.kOptimal":
        print(f"{name}: FAILED (status={status}), skipping")
        continue

    sol = h.getSolution()
    x = sol.col_value
    sum_x = sum(x)
    obj = h.getInfoValue("objective_function_value")[1]
    cut_rhs = math.floor(sum_x) + 1

    print(f"{name}: obj={obj:.4f}  sum_x={sum_x:.4f}  cut: sum_x >= {cut_rhs}")

    # Add the cut: sum(x_i for all i) >= cut_rhs
    n = h.getNumCol()
    idx = np.arange(n, dtype=np.int32)
    val = np.ones(n, dtype=np.float64)
    h.addRow(float(cut_rhs), 1e30, n, idx, val)

    # Write modified model
    h.writeModel(dst)

print("Done.")
