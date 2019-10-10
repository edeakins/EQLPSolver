#include "HDual.h"
#include "HTimer.h"
#include "HTester.h"
#include "HPrimal.h"

#include <set>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

void solvePlain(const char *filename);

int main(int argc, char **argv) {
        solvePlain(argv[1]);
    return 0;
}

void solvePlain(const char *filename) {
    // Initialize Model
    HModel model;
    model.intOption[INTOPT_PRINT_FLAG] = 1;
    model.intOption[INTOPT_SCALE_FLAG] = 0;
    // Do initial equitable partition and solve aggregated model
    model.setup(filename);
    HDual solver;
    HPrimal primalSolver;
    solver.solve(&model);
    for (int i = 0; i < model.basicIndex.size(); ++i){
	    	if (model.basicIndex[i] >= model.aggNumCol){
	    		cout << model.aggRowIdx[model.basicIndex[i] - model.aggNumCol] << " has value: " << model.baseValue[i] << endl;
	    		continue;
	    	}
    		cout << model.basicIndex[i] << " has value: " << model.baseValue[i] << endl;
	    }
    // Print the results -- writing functionality does not appear to be working
    //model.printResult();
    //model.writePivots("p");
    // Test discretizing function and resolve
    while(!model.discrete()){
	    model.isolate(model.iso);
	    model.build();
	    cout << "masterIter: " << model.masterIter << endl;
	    //model.aggregateCT();
	    //model.initCost();
	    //model.initValue();
	   	primalSolver.solvePhase2(&model);
	    for (int i = 0; i < model.basicIndex.size(); ++i){
	    	if (model.basicIndex[i] >= model.aggNumCol){
	    		cout << model.basicIndex[i] << " has value: " << model.baseValue[i] << endl;
	    		continue;
	    	}
    		cout << model.basicIndex[i] << " has value: " << model.baseValue[i] << endl;
	    }
	}

	// cout << model.oldNumCols << endl;
	// cin.get();
 //    for (int i = 0; i < model.basicVars.size(); ++i){
 //    	if (model.basicVars[i] && i < model.numCol)
 //    		cout << "var: " << i << " is basic." << endl;
 //    }
 //    cout << "\n" << endl;

    

    // for (int i = 0; i < model.basicVars.size(); ++i){
    // 	if (model.basicVars[i])
    // 		cout << "var: " << i << " is basic." << endl;
    // }
    // cout << "\n" << endl;
    
    // model.printResult();
    // for (int i = 0; i < model.aggNumCol; ++i){
    //     cout << "var color: " << model.aggColIdx[i] << endl;
    //     for (int j = model.aggAstart[i]; j < model.aggAstart[i + 1]; ++j)
    //         cout << model.aggAvalue[j] << endl;
    // }
    // for (int i  = 0; i < model.residuals.size(); ++i)
    //     cout << model.residuals[i] << endl;
    //solver.solve(&model);
    // for (int i = 0; i < model.aggNumCol; ++i){
    //     cout << "var: " << i << " obj: " << model.aggColCost[i] << endl;;
    // }
}

