\ File written by Highs .lp filereader
minimize
 obj: -1 x1 -55 x2 -330 x3 -462 x4 -165 x5 -11 x6 +0 x7 +0 x8 +0 x9 +0 x10 +0 x11 
st
 con0: +11 x0 +110 x1 = +11
 con1: +495 x1 +1320 x2 = +165
 con2: +2310 x2 +2772 x3 = +462
 con3: +2310 x3 +1320 x4 = +330
 con4: +495 x4 +110 x5 = +55
 con5: +11 x5 = +1
 con6: +1 x0 -1 x1 -1 x6 = +0
 con7: +1 x0 -1 x2 -1 x7 = +0
 con8: +1 x0 -1 x3 -1 x8 = +0
 con9: +1 x0 -1 x4 -1 x9 = +0
 con10: +1 x0 -1 x5 -1 x10 = +0
bounds
 +0 <= x0 <= +1
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
 +0 <= x4 <= +1
 +0 <= x5 <= +1
 +0 <= x6 <= +0
 +0 <= x7 <= +0
 +0 <= x8 <= +0
 +0 <= x9 <= +0
 +0 <= x10 <= +0
binary
general
semi
end
