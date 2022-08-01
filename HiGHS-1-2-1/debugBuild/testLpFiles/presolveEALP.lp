\ File written by Highs .lp filereader
min
 obj: +4 x1 +4 x2 +1 x3 
st
 con1: +8 x1 +8 x2 +4 x3 = +4
 con2: +4 x1 +1 x3 = +1
 con3: +1 x1 -1 x3 = +0
bounds
 +0 <= x1 <= +1
 +0 <= x2 <= +1
 +0 <= x3 <= +1
bin
gen
semi
end
