\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  
Subject To
 R0: 2x_2 1x_3 1x_4 1x_5 >= 1
 R1: 1x_0 1x_1 2x_2 2x_3 1x_4 2x_5 >= 2
 R2: 1x_0 1x_1 2x_3 1x_4 >= 1 
 R3: 1x_0 1x_2 2x_3 1x_5 >= 1
 R4: 1x_1 1x_2 1x_3 2x_5 >= 1
 R5: 1x_0 1x_1 2x_2 1x_4 >= 1
 R6: 1x_0 1x_1 1x_4 2x_5 >= 1
 R7: x_1 - x_3 - r_13 = 0
 R8: x_1 - x_4 - r14 = 0
Bounds
 x_0 <= 1
 x_1 <= 1
 x_2 <= 1
 x_3 <= 1
 x_4 <= 1
 x_5 <= 1
 r_13 = 0
 r_14 = 0
End
