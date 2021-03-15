import os
from gurobipy import *
import pandas as pd

timeDict = {}
wdir = "./LPGeneration/benchmarkInstances/"
for f in os.listdir(wdir):
    if f.endswith(".mps"):
        m = read(wdir + f)
        m.Params.Threads = 1
        m.optimize()
        timeDict[f] = m.Runtime

timeDF = pd.DataFrame.from_dict(timeDict, orient="index")
timeDF.to_csv("gurobi_times.csv", header="False")