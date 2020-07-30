\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  
Subject To
R0: 0.5 x_5 + 0.5 x_8 = 0
R1: x_0 = 0.5
R2: x_1 = 0.5
R3: x_0 + x_2 + x_3 + x_5 + x_6 = 1
R4: x_0 + x_1 + x_2 + x_4 + x_7 >= 1
R5: 0.5 x_2 + 0.5 x_7 = 0.5
R6: x_3 - x_6 = 0
R7: x_5 - x_8 = 0
R8: x_2 - x_7 = 0
R9: x_1 - x_2 - x_4 + x_5 = 0
Bounds
 x_0 <= 1
 x_1 <= 1
 x_2 <= 1
 x_3 <= 1
 x_4 <= 1
 x_5 <= 1
 x_6 <= 1
 x_7 <= 1
 x_8 <= 1
End