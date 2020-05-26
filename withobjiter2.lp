\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  35 C0 + C1 + C2 + 21 C3 + 21 C4 + 105 C5 + 35 C6 + 21 C7 + 63 C8 + 7 C9
   + 3 C10 + 7 C11 + C12 + 21 C13 + 35 C14 + 7 C15
Subject To
 R0: 5 C0 + 3 C8 + 2 C11 = 1
 R1: C1 + 6 C3 + 3 C4 = 1
 R2: C2 + 3 C4 + 6 C7 = 1
 R3: C1 + C2 + 7 C4 + C12 = 1
 R4: C3 + 2 C4 + 5 C5 + C7 + C13 = 1
 R5: 3 C5 + 4 C6 + 3 C7 = 1
 R6: 4 C0 + 3 C3 + 3 C5 = 1
 R7: C0 + 4 C5 + C6 + 3 C8 + C14 = 1
 R8: 5 C6 + 3 C8 + 2 C9 = 1
 R9: 6 C8 + C9 + C10 + C11 + C15 = 1
 R10: 7 C9 + 3 C10 = 1
 R11: 3 C10 + 7 C11 >= 1
 R12: 3 C4 + C12 + 6 C13 = 1
 R13: 3 C5 + 3 C13 + 4 C14 = 1
 R14: 3 C8 + 5 C14 + 2 C15 = 1
 R15: 3 C10 + 7 C15 = 1
Bounds
 C0 <= 1
 C1 <= 1
 C2 <= 1
 C3 <= 1
 C4 <= 1
 C5 <= 1
 C6 <= 1
 C7 <= 1
 C8 <= 1
 C9 <= 1
 C10 = 0
 C11 <= 1
 C12 <= 1
 C13 <= 1
 C14 <= 1
 C15 <= 1
End
