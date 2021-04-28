\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9
Subject To
 R0: x_1 + x_2 + x_3 + x_4 + x_7 >= 1
 R1: x_1 + x_2 + x_3 + x_5 + x_8 >= 1
 R2: x_1 + x_2 + x_3 + x_6 + x_9 >= 1
 R3: x_1 + x_4 + x_5 + x_6 + x_7 >= 1
 R4: x_2 + x_4 + x_5 + x_6 + x_8 >= 1
 R5: x_3 + x_4 + x_5 + x_6 + x_9 >= 1
 R6: x_1 + x_4 + x_7 + x_8 + x_9 >= 1
 R7: x_2 + x_5 + x_7 + x_8 + x_9 >= 1
 R8: x_3 + x_6 + x_7 + x_8 + x_9 >= 1
 R9: x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 = 2
 R10: +1 x_5 -1 x_7 -1 x_10 = +0
 R11: +1 x_6 -1 x_8 -1 x_11 = +0
 R12: +1 x_6 -1 x_9 -1 x_12 = +0
Bounds
 x_1 <= 1
 x_2 <= 1
 x_3 <= 1
 x_4 <= 1
 x_5 <= 1
 x_6 <= 1
 x_7 <= 1
 x_8 <= 1
 x_9 <= 1
 x_10 <= 0
 x_11 <= 0
 x_12 <= 0
End
