#include "HPrimal.h"
#include "HConst.h"

#include <cstdio>
#include <iostream>
using namespace std;

void HPrimal::solvePhase2(HModel *ptr_model) {
    // cout << "numRegCol: " << model->oldNumCols << endl;
    // cin.get();
    // Copy size
    model = ptr_model;
    numCol = model->getNumCol();
    numRow = model->getNumRow();
    numTot = model->getNumTot();
    model->computeFactor();

    // Setup update limits
    limitUpdate = min(100 + numRow / 100, 1000);
    countUpdate = 0;

    // Setup local vectors
    column.setup(numRow);
    row_ep.setup(numRow);
    row_ap.setup(numCol);
    columnDensity = 0;
    row_epDensity = 0;

    // for (int i = 0; i < model->aggNumRow; ++i){
    //     cout << "Row: " << i << " lower: " << model->aggRowLower[i] << " upper: " << model->aggRowUpper[i] << endl;
    // }
    // for (int i = model->oldNumCols; i < model->aggNumCol; ++i){
    //     cout << "var: " << i << " lower: " << model->aggColLower[i] << " upper: " << model->aggColUpper[i] << endl;
    // }
    // cin.get();

   	// Unfolding - Deakins, Ostrowski, Kneuven, Schrock
    if (model->masterIter - 1){
    	primalRebuild();
        model->printMessage("unfolding-start");
        primalChooseLinker();
        primalCollect();
        return;
    }

    // Setup other buffers - Huha
    model->printMessage("primal-start");
    for (;;) {
        primalRebuild();
        cin.get();
        for (;;) {
            primalChooseColumn();
            if (columnIn == -1)
                break;

            primalChooseRow();
            if (rowOut == -1) {
                break;
            }
            primalUpdate();
            if (invertHint)
                break;
        }

        // Fresh factorization required
        if (countUpdate == 0)
            break;
    }
    if (columnIn == -1) {
        model->printMessage("primal-optimal");
        model->printMessage("problem-optimal");
        model->reportStatus(LP_Status_Optimal);
    } else {
        model->printMessage("primal-unbounded");
        model->reportStatus(LP_Status_Unbounded);
    }
}

// Choose Linker
void HPrimal::primalChooseLinker(){
	 for (int i = 0; i < model->residuals.size(); ++i){
    	// if (model->basicResiduals[i]){
     //        // cout << i << endl;
    	// 	continue;
     //    }
        columnIn = model->residuals[i]; 
        //cout << "columnIN: " << columnIn << endl;
        cout << "columnIn: " << columnIn << endl;   
        //model->basicVars[columnIn] = true;
        model->workCost[columnIn] = 1;
        model->workUpper[columnIn] = +HSOL_CONST_INF;
        model->workLower[columnIn] = -HSOL_CONST_INF; 
        model->workValue[columnIn] = 0;
        model->nonbasicMove[columnIn] = 1;
        primalChooseRow();
        primalUpdate();
        primalRebuild();
        //setNow();
        // for (int j = 0; j < model->aggRowIdx.size(); ++j){
        //     cout << "row: " << j << " = " << model->aggRowUpper[j] << endl;
        // }
        // cin.get();
        model->workCost[columnIn] = 0;
    }
}

void HPrimal::setNow(){
    vector<bool> activeCon(model->oldNumRows, true);
    vector<bool> activeVar(model->oldNumCols, true);
    int slackIdx = 0;
    for (int i = 0; i < model->basicIndex.size(); ++i){
        if (model->basicIndex[i] < model->aggNumCol){
            if (model->basicIndex[i] < model->oldNumCols){
                //cout << "i: " << i << endl;
                if (fabs(model->baseValue[i]) < 1e-5){
                    continue;
                }
                else{
                    activeVar[i] = false;
                }
                //model->boundedVariables[model->aggColIdx[model->basicIndex[i]]] = upperBound + lowerBound;
            }
        }
        else{
            slackIdx = model->basicIndex[i] - model->aggNumCol;
            activeCon[slackIdx] = false;
            //model->activeConstraints[slackIdx] = false;
        }
    }
    for (int i = 0; i < activeCon.size(); ++i){
        if (activeCon[i]){
        model->aggRowUpper[i] = min(model->aggRowUpper[i], model->aggRowLower[i]);
        model->aggRowLower[i] = min(model->aggRowUpper[i], model->aggRowLower[i]);
    }
    }
}

void HPrimal::primalCollect(){
    vector<bool> active(model->aggNumRow, true);
    vector<bool> nonBasic(model->oldNumCols, true);
	int slackIdx = 0;
    for (int i = 0; i < model->basicIndex.size(); ++i){
        //cout << "var: " << model->basicIndex[i] << " = " << model->baseValue[i] << endl;
        if (model->basicIndex[i] < model->aggNumCol){
        	if (model->basicIndex[i] < model->oldNumCols){
                if (fabs(model->baseValue[i]) < 1e-5){
                    //cout << "degenerate: " << i << endl;
                    continue;
                }
                nonBasic[model->basicIndex[i]] = false;
                model->prevBasicColor[model->aggColIdx[model->basicIndex[i]]] = true;
                model->prevBasicValue[model->aggColIdx[model->basicIndex[i]]] = model->baseValue[i];
                bool upperBound = model->baseValue[i] == model->aggColUpper[model->basicIndex[i]];
                bool lowerBound = model->baseValue[i] == model->aggColUpper[model->basicIndex[i]];
                //model->boundedVariables[model->aggColIdx[model->basicIndex[i]]] = upperBound + lowerBound;
            }
        }
        // else{
        //     slackIdx = model->basicIndex[i] - model->aggNumCol;
        //     cout << slackIdx << endl;
        //     cout << "row: " << model->aggRowIdx[slackIdx] << endl;
        //     cin.get();
        //     model->prevBasicColor[model->aggRowIdx[slackIdx]] = true;
        //     model->prevBasicValue[model->aggRowIdx[slackIdx]] = model->baseValue[i];
        //     //model->activeConstraints[slackIdx] = false;
        // }
   	}
    // for (int i = 0; i < nonBasic.size(); ++i){
    //     if (nonBasic[i])
    //         model->prevBasicValue[i] = model->workValue[i];
    // }
    model->collectActiveConstraints();
    for (int i = 0; i < model->oldNumRows; ++i){
        if (!model->activeConstraints[i]){
            model->prevBasicColor[model->aggRowIdx[i]] = true;
            model->prevBasicValue[model->aggRowIdx[i]] = model->baseValue[i];
        }
    }

    //cin.get();

   	// for (int i = 0; i < model->activeConstraints.size(); ++i)
   	// 	cout << model->activeConstraints[i] << endl;
   	// cin.get();
}

void HPrimal::primalRebuild() {
    model->reportPivots(-1, -1, 0); // Indicate REINVERT

    // Rebuild model->factor - only if we got updates
    invertHint = 0;
    if (countUpdate > 0) {
        model->computeFactor();
        countUpdate = 0;
    }

    model->computeDual();
    model->computePrimal();
    model->computeObject();
    model->printObject();
}

void HPrimal::primalChooseColumn() {
    columnIn = -1;
    double bestInfeas = 0;
    const int *jFlag = model->getNonbasicFlag();
    const int *jMove = model->getNonbasicMove();
    double *workDual = model->getWorkDual();
    const double *workLower = model->getWorkLower();
    const double *workUpper = model->getWorkUpper();
    const double dualTolerance = model->dblOption[DBLOPT_DUAL_TOL];

    for (int iCol = 0; iCol < numTot; iCol++) {
        if (jFlag[iCol] && fabs(workDual[iCol]) > dualTolerance) {
            // Always take free
            // TODO: if we found free,
            // Then deal with it in dual phase 1
            if (workLower[iCol] == -HSOL_CONST_INF
                    && workUpper[iCol] == HSOL_CONST_INF) {
                columnIn = iCol;
                break;
            }
            // Then look at dual infeasible
            if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
                if (bestInfeas < fabs(workDual[iCol])) {
                    bestInfeas = fabs(workDual[iCol]);
                    columnIn = iCol;
                }
            }
        }
    }
}

void HPrimal::primalChooseRow() {
    const double *baseLower = model->getBaseLower();
    const double *baseUpper = model->getBaseUpper();
    double *baseValue = model->getBaseValue();
    const double primalTolerance = model->dblOption[DBLOPT_PRIMAL_TOL];

    // Compute pivot column
    column.clear();
    column.packFlag = true;
    model->getMatrix()->collect_aj(column, columnIn, 1);
    model->getFactor()->ftran(column, columnDensity);
    columnDensity = 0.95 * columnDensity + 0.05 * column.count / numRow;

    // Initialize
    rowOut = -1;

    // Choose column pass 1
    double alphaTol = countUpdate < 10 ? 1e-9 : countUpdate < 20 ? 1e-8 : 1e-7;
    int *jMove = model->getNonbasicMove();
    //jMove[columnIn] = 1;
    int moveIn = jMove[columnIn];
    if (moveIn == 0) {
        // If there's still free in the N
        // We would report not-solved
        // Need to handle free
    }
    double relaxTheta = 1e100;
    for (int i = 0; i < column.count; i++) {
        int index = column.index[i];
        double alpha = column.array[index] * moveIn;
        if (alpha > alphaTol) {
            double relaxSpace = baseValue[index] - baseLower[index]
                    + primalTolerance;
            if (relaxSpace < relaxTheta * alpha)
                relaxTheta = relaxSpace / alpha;
        } else if (alpha < -alphaTol) {
            double relaxSpace = baseValue[index] - baseUpper[index]
                    - primalTolerance;
            if (relaxSpace > relaxTheta * alpha)
                relaxTheta = relaxSpace / alpha;
        }
    }

    // Choose column pass 2
    double bestAlpha = 0;
    for (int i = 0; i < column.count; i++) {
        int index = column.index[i];
        double alpha = column.array[index] * moveIn;
        if (alpha > alphaTol) {
            double tightSpace = baseValue[index] - baseLower[index];
            if (tightSpace < relaxTheta * alpha) {
                if (bestAlpha < alpha) {
                    bestAlpha = alpha;
                    rowOut = index;
                }
            }
        } else if (alpha < -alphaTol) {
            double tightSpace = baseValue[index] - baseUpper[index];
            if (tightSpace > relaxTheta * alpha) {
                if (bestAlpha < -alpha) {
                    bestAlpha = -alpha;
                    rowOut = index;
                }
            }
        }
    }
    // cout << "rowOut: " << rowOut << endl;
}

void HPrimal::primalUpdate() {
    int *jMove = model->getNonbasicMove();
   // jMove[columnIn] = 1;
    double *workDual = model->getWorkDual();
    const double *workLower = model->getWorkLower();
    const double *workUpper = model->getWorkUpper();
    const double *baseLower = model->getBaseLower();
    const double *baseUpper = model->getBaseUpper();
    double *workValue = model->getWorkValue();
    double *baseValue = model->getBaseValue();
    const double primalTolerance = model->dblOption[DBLOPT_PRIMAL_TOL];

    // Compute thetaPrimal
    int moveIn = jMove[columnIn];
    int columnOut = model->getBaseIndex()[rowOut];
    double alpha = column.array[rowOut];
    double thetaPrimal = 0;
    if (alpha * moveIn > 0) {
        // Lower bound
        thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
    } else {
        // Upper bound
        thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
    }

    // 1. Make sure it is inside bounds or just flip bound
    double lowerIn = workLower[columnIn];
    double upperIn = workUpper[columnIn];
    double valueIn = workValue[columnIn] + thetaPrimal;
    bool flipped = false;
    if (jMove[columnIn] == 1) {
        if (valueIn > upperIn + primalTolerance) {
            // Flip to upper
            workValue[columnIn] = upperIn;
            thetaPrimal = upperIn - lowerIn;
            flipped = true;
            jMove[columnIn] = -1;
        }
    } else if (jMove[columnIn] == -1) {
        if (valueIn < lowerIn - primalTolerance) {
            // Flip to lower
            workValue[columnIn] = lowerIn;
            thetaPrimal = lowerIn - upperIn;
            flipped = true;
            jMove[columnIn] = 1;
        }
    }

    for (int i = 0; i < column.count; i++) {
        int index = column.index[i];
        baseValue[index] -= thetaPrimal * column.array[index];
    }

    // If flipped, then no need touch the pivots
    if (flipped) {
        return;
    }

    // Pivot in
    int sourceOut = alpha * moveIn > 0 ? -1 : 1;
    model->updatePivots(columnIn, rowOut, sourceOut);

    baseValue[rowOut] = valueIn;

    // Check for any possible infeasible
    for (int iRow = 0; iRow < numRow; iRow++) {
        if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
            invertHint = 1;
        } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
            invertHint = 1;
        }
    }

    // 2. Now we can update the dual
    row_ep.clear();
    row_ap.clear();
    row_ep.count = 1;
    row_ep.index[0] = rowOut;
    row_ep.array[rowOut] = 1;
    row_ep.packFlag = true;
    model->getFactor()->btran(row_ep, row_epDensity);
    model->getMatrix()->price_by_row(row_ap, row_ep);
    row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / numRow;

    double thetaDual = workDual[columnIn] / alpha;
    for (int i = 0; i < row_ap.count; i++) {
        int iCol = row_ap.index[i];
        workDual[iCol] -= thetaDual * row_ap.array[iCol];
    }
    for (int i = 0; i < row_ep.count; i++) {
        int iGet = row_ep.index[i];
        int iCol = iGet + numCol;
        workDual[iCol] -= thetaDual * row_ep.array[iGet];
    }

    // Dual for the pivot
    workDual[columnIn] = 0;
    workDual[columnOut] = -thetaDual;

    // Update model->factor basis
    model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
    model->updateMatrix(columnIn, columnOut);
    model->collectActiveConstraints();
    // model->updateActiveConstraints();

    if (++countUpdate >= limitUpdate)
        invertHint = true;
    model->reportPivots(columnIn, columnOut, alpha);
    // for (int i = 0; i < column.array.size(); ++i)
    // 	cout << "idx: " << column.index[i] << " = " << column.array[i] << endl;
    cout << "columnOut: " << columnOut << endl;
    cin.get();
}
