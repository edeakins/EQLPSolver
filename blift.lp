 \ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
 r01
Subject To
R0: 1/9 x0 + 4/9 x1 + 4/9 x2 = 2/9
R1: 4/9 x1 + 4/9 x2 - 1/9 s0 - 4/9 s1 - 4/9 s2 = 1/9
R2: s0 - s1 >= 0
R3: s0 - s2 >= 0
R4: x0 - x1 - r01 = 0
R5: x0 - x2 - r02 = 0
Bounds
x0 <= 1
x1 <= 1
x2 <= 1
r02 = 0
s0 free
s1 free
s2 free
End