\ File written by Highs .lp filereader
minimize
 obj: -1 x1 -576 x2 -14 x3 -84 x4 -280 x5 -560 x6 -672 x7 +0 x8 +0 x9 +0 x10 +0 x11 +0 x12 +0 x13 
st
 con1: +14 x1 +28 x3 +168 x4 = +14
 con2: +5184 x2 +2688 x7 = +576
 con3: +1 x1 +14 x3 = +1
 con4: +168 x3 +252 x4 +840 x5 = +84
 con6: +3456 x2 +3360 x6 +4032 x7 = +672
 con7: +840 x4 +1120 x5 +2240 x6 = +280
bounds
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
 +0 <= x4 <= +1
 +0 <= x5 <= +1
 +0 <= x6 <= +1
 +0 <= x7 <= +1
 +0 <= x8 <= +0
 +0 <= x9 <= +0
 +0 <= x10 <= +0
 +0 <= x11 <= +0
 +0 <= x12 <= +0
 +0 <= x13 <= +0
binary
general
semi
end