\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  r_14
Subject To
 R0: 2 x_2 + 1 x_3 + 1 x_4 + 1 x_5 = 1
 R1: x_0 = .5
 R2: 1 x_0 + 1 x_1 + 2 x_3 + 1 x_4 = 1 
 R3: 1 x_0 + 1 x_2 + 2 x_3 + 1 x_5 >= 1
 R4: 1 x_0 + 1 x_1 + 2 x_2 + 1 x_4 >= 1
 R5: 1 x_0 + 1 x_1 + 1 x_4 + 2 x_5 >= 1
 R6: x_1 - x_2 - x_4 + x_5 = 0
 R7: x_1 - x_3 - r_13 = 0
 R8: x_1 - x_4 - r_14 = 0
 R9: x_3 = 0
Bounds
 x_0 <= 1
 x_1 <= 1
 x_2 <= 1
 x_3 <= 1
 x_4 <= 1
 x_5 <= 1
End