NAME
ROWS
 N  COST
 E  R0      
 E  R1      
 E  R2      
 G  R3      
 G  R4      
 G  R5      
 E  R6      
 E  R7      
 E  R8      
COLUMNS
    C0        COST      1
    C0        R1        1
    C0        R3        1
    C0        R4        1
    C0        R5        1
    C1        COST      1
    C1        R4        1
    C1        R5        1
    C1        R7        1
    C1        R8        1
    C2        COST      2
    C2        R3        1
    C2        R4        2
    C3        COST      2
    C3        R3        2
    C3        R7        -1
    C4        COST      1
    C4        R4        1
    C4        R5        1
    C4        R8        -1
    C5        COST      2
    C5        R3        1
    C5        R5        2
    C6        R7        -1
    C7        R8        -1
RHS
    RHS_V     R0        0.25
    RHS_V     R1        0.5
    RHS_V     R2        0.375
    RHS_V     R3        1
    RHS_V     R4        1
    RHS_V     R5        1
BOUNDS
 UP BOUND     C0        1
 UP BOUND     C1        1
 UP BOUND     C2        1
 UP BOUND     C3        1
 UP BOUND     C4        1
 UP BOUND     C5        1
 FX BOUND     C6        0
 FX BOUND     C7        0
ENDATA
