\ Fractional matching LP of K4 (complete graph on 4 vertices), x_i indexing.
\ x1=(1,2) x2=(1,3) x3=(1,4) x4=(2,3) x5=(2,4) x6=(3,4)
\ Constraint c_v = sum of edges incident to vertex v <= 1.
\ 6 variables > 4 constraints. Each x_i is in exactly 2 constraints; each
\ constraint contains exactly 3 variables -> biregular, uniform data.
Maximize
 obj: x1 + x2 + x3 + x4 + x5 + x6
Subject To
 c1: x1 + x2 + x3 <= 1
 c2: x1 + x4 + x5 <= 1
 c3: x2 + x4 + x6 <= 1
 c4: x3 + x5 + x6 <= 1
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
End
