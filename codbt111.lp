\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0 + x1 + x2 + x3 + x4 + x5
Subject To
 R0: x0 + x1 + x2 + x3 >= 1
 R1: x0 + x1 + x2 + x4 >= 1
 R2: x0 + x1 + x2 + x5 >= 1
 R3: x0 + x3 + x4 + x5 >= 1
 R4: x1 + x3 + x4 + x5 >= 1
 R5: x2 + x3 + x4 + x5 >= 1
Bounds
 x0 <= 1
 x1 <= 1
 x2 <= 1
 x3 <= 1
 x4 <= 1
 x5 <= 1
End
