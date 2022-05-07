\ File written by Highs .lp filereader
min
 obj: +4 x1 +4 x2 +1 x3 +0 x4 +0 x5 
st
 con1: +8 x1 +8 x2 +4 x3 = +4
 con2: +8 x1 +12 x2 = +4
 con3: +4 x1 +1 x3 = +1
 con4: +1 x1 -1 x2 -1 x4 = +0
 con5: +1 x1 -1 x3 -1 x5 = +0
bounds
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
 x4 free
 x5 free
end
