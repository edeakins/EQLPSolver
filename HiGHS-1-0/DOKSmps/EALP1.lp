\ File written by Highs .lp filereader
minimize
 obj: +1 x1 +25 x2 +100 x3 +100 x4 +25 x5 +1 x6 +0 x7 +0 x8 +0 x9 +0 x10 +0 x11 
st
 con0: +10 x0 +150 x1 +300 x2 +100 x3 = +10
 con1: +100 x1 +600 x2 +600 x3 +100 x4 = +25
 con2: +100 x2 +300 x3 +150 x4 +10 x5 = +10
 con3: +1 x0 -1 x1 -1 x6 = +0
 con4: +1 x0 -1 x2 -1 x7 = +0
 con5: +1 x0 -1 x3 -1 x8 = +0
 con6: +1 x0 -1 x4 -1 x9 = +0
 con7: +1 x0 -1 x5 -1 x10 = +0
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
