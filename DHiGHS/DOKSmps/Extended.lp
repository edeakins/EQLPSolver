\ File written by Highs .lp filereader
minimize
 obj: +1 x1 +1 x2 +1 x3 +0 x4 +0 x5 
st
 con1: +1 x1 +1 x2 = +1
 con2: +1 x2 +1 x3 = +1
 con3: +1 x1 +1 x3 = +1
 con4lo: +1 x1 +1 x2 +1 x3 >= +1
 con5: +1 x1 -1 x2 -1 x4 = +0
 con6: +1 x1 -1 x3 -1 x5 = +0
bounds
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
 +0 <= x4 <= +0
 +0 <= x5 <= +0
binary
general
semi
end
