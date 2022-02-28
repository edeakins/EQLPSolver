obj: min 
    - y0 - y1 - y2 - y3 - y4 - y5 - y6 - y7 - y8 + y9 + y10 + y11 + y12 + y13 + y14 + y15 + y16 + y17
s.t.
row0: y0 + y1 + y2 + y3           + y6           - y9                                                 + z0                                         = 1
row1: y0 + y1 + y2      + y4           + y7           - y10                                                + z1                                    = 1
row2: y0 + y1 + y2           + y5           + y8            - y11                                               + z2                               = 1
row3: y0           + y3 + y4 + y5 + y6                            - y12                                              + z3                          = 1
row4:      y1      + y3 + y4 + y5      + y7                             - y13                                             + z4                     = 1
row5:           y2 + y3 + y4 + y5           + y8                              - y14                                            + z5                = 1
row6: y0           + y3           + y6 + y7 + y8                                    - y15                                           + z6           = 1
row7:      y1           + y4      + y6 + y7 + y8                                          - y16                                          + z7      = 1
row8:           y2           + y5 + y6 + y7 + y8                                                - y17                                         + z8 = 1

Bounds:
yi, zi <= inf
