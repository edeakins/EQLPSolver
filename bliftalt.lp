\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0 + 4 x1 + 4 x2
Subject To
 R0: x0 + 4 x1 + 4 x2 = 2
 R1: - s0 = 1
 R2: s0 - s1 = 0
 R3: s0 - s2 = 0
 R4: x0 - x1 - r01 = 0
 R5: x0 - x2 - r02 = 0
Bounds
 x0 <= 1
 x1 <= 1
 x2 <= 1
 s0 free
 s1 free
 s2 free
 r01 = 0
 r02 = 0
End
