#include "HDual.h"
#include "HTimer.h"
#include "HTester.h"

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
    model.setup(filename);
    // Testing
    
    // Set solver and solve the model
    HDual solver;
    solver.solve(&model);
    // Print the results -- writing functionality does not appear to be working
    model.printResult();
    model.writePivots("p");
    model.discrete();
    model.isolate(model.iso);
    
    model.build();
    solver.solve(&model);
    // for (int i = 0; i < model.aggNumCol; ++i){
    //     cout << "var: " << i << " obj: " << model.aggColCost[i] << endl;;
    // }
}

