\ File written by Highs .lp filereader
minimize
 obj: +1 x1 +25 x2 +100 x3 +100 x4 +25 x5 +1 x6 +0 x7 +0 x8 +0 x9 +0 x10 +0 x11 
st
 con1: +10 x1 +150 x2 +300 x3 +100 x4 = +10
 con2: +100 x2 +600 x3 +600 x4 +100 x5 = +25
 con3: +100 x3 +300 x4 +150 x5 +10 x6 = +10
 con4: +1 x1 -1 x2 -1 x7 = +0
 con5: +1 x1 -1 x3 -1 x8 = +0
 con6: +1 x1 -1 x4 -1 x9 = +0
 con7: +1 x1 -1 x5 -1 x10 = +0
 con8: +1 x1 -1 x6 -1 x11 = +0
bounds
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
 +0 <= x4 <= +1
 +0 <= x5 <= +1
 +0 <= x6 <= +1
 +0 <= x7 <= +0
 +0 <= x8 <= +0
 +0 <= x9 <= +0
 +0 <= x10 <= +0
 +0 <= x11 <= +0
binary
general
semi
end
