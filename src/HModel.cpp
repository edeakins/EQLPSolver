#include "HModel.h"
#include "HConst.h"
#include "HTimer.h"

#include <cctype>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
#include <set>
#include <list>
#include <tuple>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <numeric>
using namespace std;

HModel::HModel() {
    intOption[INTOPT_PRINT_FLAG] = 0; // no print
    intOption[INTOPT_TRANSPOSE_FLAG] = 0; // no transpose
    intOption[INTOPT_SCALE_FLAG] = 1; // do scale
    intOption[INTOPT_TIGHT_FLAG] = 1; // do tight
    intOption[INTOPT_PERMUTE_FLAG] = 0; // no permute
    intOption[INTOPT_PERTURB_FLAG] = 1; // do perturb

    dblOption[DBLOPT_PRIMAL_TOL] = 1e-7;
    dblOption[DBLOPT_DUAL_TOL] = 1e-7;
    dblOption[DBLOPT_PAMI_CUTOFF] = 0.95;

    strOption[STROPT_PARTITION_FILE] = "";
}

void HModel::setup(const char *filename) {
    // Setup timer
    timer.reset();
    totalTime = 0;
    problemStatus = -1;
    numberIteration = 0;
    masterIter = 0;
    modelName = filename;

    // Load the model
    setup_loadMPS(filename);

    // Check size
    numTot = numCol + numRow;
    if (numRow == 0)
        return;
    
 	initEQs();
    // init and build all working pieces of model
    build();
}
// Split up setup func so that we could re initialize agg model attributes
void HModel::build(){
	startingBasis.clear();
	// Compute EP refinement
    equitable();
    // Update Master iter count
    masterIter++;
    // Setup random buffers: shuffle the break
    intBreak.resize(aggNumTot);
    for (int i = 0; i < aggNumTot; i++)
        intBreak[i] = i;
    for (int i = aggNumTot - 1; i >= 1; i--) {
        int j = random.intRandom() % (i + 1);
        swap(intBreak[i], intBreak[j]);
    }
    dblXpert.resize(aggNumTot);
    for (int i = 0; i < aggNumTot; i++)
        dblXpert[i] = random.dblRandom();

    // Setup other part
    setup_transposeLP();
    setup_scaleMatrix();
    //    setup_tightenBound();
    //    setup_shuffleColumn();
    setup_allocWorking();

    // Initialize the values
    initCost();
    initBound();
    initValue();

    // Save the input time
    totalTime += timer.getTime();
}

bool load_mpsLine(FILE *file, int lmax, char *line, char *flag, double *data) {
    int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;

    // check the buffer
    if (flag[1]) {
        flag[1] = 0;
        memcpy(&data[2], &line[F5], 8);
        data[0] = atof(&line[F6]);
        return true;
    }

    // try to read some to the line
    for (;;) {
        // Line input
        fgets(line, lmax, file);

        // Line trim   -- to delete tailing white spaces
        int lcnt = strlen(line) - 1;
        while (isspace(line[lcnt]) && lcnt >= 0)
            lcnt--;
        if (lcnt <= 0 || line[0] == '*')
            continue;

        // Line expand -- to get data easier
        lcnt++;
        while (lcnt < F4)
            line[lcnt++] = ' '; // For row and bound row name
        if (lcnt == F4)
            line[lcnt++] = '0'; // For bound value
        line[lcnt] = '\0';

        // Done with section symbol
        if (line[0] != ' ') {
            flag[0] = line[0];
            return false;
        }

        // Read major symbol & name
        flag[0] = line[F1 + 1] == ' ' ? line[F1] : line[F1 + 1];
        memcpy(&data[1], &line[F2], 8);

        // Read 1st minor name & value to output
        memcpy(&data[2], &line[F3], 8);
        data[0] = atof(&line[F4]);

        // Keep 2nd minor name & value for future
        if (lcnt > F5)
            flag[1] = 1;
        break;
    }

    return true;
}

void HModel::setup_loadMPS(const char *filename) {
    // MPS file buffer
    numRow = 0;
    numCol = 0;
    objOffset = 0;
    FILE *file = fopen(filename, "r");
    if (file == 0) {
        return;
    }

    // Input buffer
    const int lmax = 128;
    char line[lmax];
    char flag[2] = { 0, 0 };
    double data[3];

    // Load NAME and ROWS
    load_mpsLine(file, lmax, line, flag, data);
    load_mpsLine(file, lmax, line, flag, data);

    vector<char> rowType;
    map<double, int> rowIndex;
    double objName = 0;
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (flag[0] == 'N' && objName == 0) {
            objName = data[1];
        } else {
            rowType.push_back(flag[0]);
            rowIndex[data[1]] = numRow++;
        }
    }

    // Load COLUMNS
    map<double, int> colIndex;
    double lastName = 0;
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (lastName != data[1]) { // New column
            lastName = data[1];
            colIndex[data[1]] = numCol++;
            colCost.push_back(0);
            Astart.push_back(Aindex.size());
        }
        if (data[2] == objName) // Cost
            colCost.back() = data[0];
        else if (data[0] != 0) {
            Aindex.push_back(rowIndex[data[2]]);
            Avalue.push_back(data[0]);
        }
    }
    Astart.push_back(Aindex.size());

    // Load model
    vector<double> RHS(numRow, 0);
    while (load_mpsLine(file, lmax, line, flag, data)) {
        if (data[2] != objName) {
            int iRow = rowIndex[data[2]];
            RHS[iRow] = data[0];
        } else {
            objOffset = data[0]; // Objective offset
        }
    }
    Rhs = RHS;

    // Load RANGES
    rowLower.resize(numRow);
    rowUpper.resize(numRow);
    if (flag[0] == 'R') {
        while (load_mpsLine(file, lmax, line, flag, data)) {
            int iRow = rowIndex[data[2]];
            if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0)) {
                rowLower[iRow] = RHS[iRow] - fabs(data[0]);
                rowUpper[iRow] = RHS[iRow];
            } else {
                rowUpper[iRow] = RHS[iRow] + fabs(data[0]);
                rowLower[iRow] = RHS[iRow];
            }
            rowType[iRow] = 'X';
        }
    }

    // Setup bounds for row without 'RANGE'
    for (int iRow = 0; iRow < numRow; iRow++) {
        switch (rowType[iRow]) {
        case 'L':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'G':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'E':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'N':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'X':
            break;
        }
    }

    // Load BOUNDS
    colLower.assign(numCol, 0);
    colUpper.assign(numCol, HSOL_CONST_INF);
    if (flag[0] == 'B') {
        while (load_mpsLine(file, lmax, line, flag, data)) {
            int iCol = colIndex[data[2]];

            switch (flag[0]) {
            case 'O': /*LO*/
                colLower[iCol] = data[0];
                break;
            case 'I': /*MI*/
                colLower[iCol] = -HSOL_CONST_INF;
                break;
            case 'L': /*PL*/
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'X': /*FX*/
                colLower[iCol] = data[0];
                colUpper[iCol] = data[0];
                break;
            case 'R': /*FR*/
                colLower[iCol] = -HSOL_CONST_INF;
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'P': /*UP*/
                colUpper[iCol] = data[0];
                if (colLower[iCol] == 0 && data[0] < 0)
                    colLower[iCol] = -HSOL_CONST_INF;
                break;
            }
        }
    }

    // Load ENDATA and close file
    fclose(file);
}

/* Deakins - TO DO: Speed this SHIT up */
// Initialize storage containers that are accesible 
// by the model class
void HModel::initStorage(){
	adjList.assign(numRow + numCol, NULL);
	C.assign(numRow + numCol, NULL);
	for (int i = 0; i < numRow + numCol; ++i){
		A.push_back(new list<int>);
	}
	SCheck.assign(numRow + numCol, false);
	mincdeg.assign(numRow + numCol, 0);
	maxcdeg.assign(numRow + numCol, 0);
	cdeg.assign(numRow + numCol, 0);
	isAdj.assign(numRow + numCol, 0);
	isolates.assign(numRow + numCol, false);
	prevBasicColor.assign(numCol + numRow, false);
	prevBasicValue.assign(numCol + numRow, false);
}

// Append to the end of linked list
void HModel::append(Node **headRef, int newData){
	Node *newNode = new Node();
	Node *last = (*headRef);
	newNode->data = newData;
	newNode->next = NULL;
	if ((*headRef) == NULL){
		newNode->prev = NULL;
		(*headRef) = newNode;
		return;
	}
	while(last->next != NULL){
		last = last->next;
	}
	last->next = newNode;
	newNode->prev = last;
}

// Append func for adjacency list of graph
void HModel::appendAdj(adjNode **headRef, int newData, double newEntry){
	adjNode *newNode = new adjNode();
	adjNode *last = (*headRef);
	newNode->label = newData;
	newNode->w = newEntry;
	newNode->next = NULL;
	if ((*headRef) == NULL){
		(*headRef) = newNode;
		return;
	}
	while(last->next != NULL){
		last = last->next;
	}
	last->next = newNode;
}

// Delete from doubly linked list function
void HModel::deleteNode(Node **headRef, Node *delet){
	if (*headRef == NULL || delet == NULL){
		return;
	}
	if (*headRef == delet){
		*headRef = delet->next;
	}
	if(delet->next != NULL){
		delet->next->prev = delet->prev;
	}
	if (delet->prev != NULL){
		delet->prev->next = delet->next;
	}
	free(delet);
	return;
}

// Delete from doubly linked list function
HModel::Node *HModel::findNode(Node **headRef, int label){
	Node *delet = (*headRef);
	while(delet != NULL){
		if (delet->data == label){
			return delet;
		}
		delet = delet->next;
	}
	return NULL;
}

// Doubly linked list searcher function
bool HModel::exists(Node **headRef, int data){
	Node *curr = (*headRef);
	while(curr != NULL){
		if (curr->data == data){
			return true;
		}
		curr = curr->next;
	}
	return false;
}

int HModel::listSize(Node **headRef){
	Node *node = (*headRef);
	int siz = 0;
	while (node != NULL){
		siz++;
		node = node->next;
	}
	return siz;
}

// Printer
void HModel::printList(Node **headRef)  
{  
	Node *node = (*headRef);
	if (!node){
		return;
	} 
    //cout << "\ncolor class: " << node->data << endl;  
    while (node != NULL)  
    {  
        cout<<" "<< node->data <<" ";  
        //last = node;  
        node = node->next;  
    }  
}

// Initialize LP as adjacency list for ease of use
void HModel::lp2Graph(){
	vector<double>AvalueCopy(Avalue.size(), 0);
	vector<double>::iterator min = min_element(Avalue.begin(), Avalue.end());
	if (*min < 0){
		for (int i = 0; i < Avalue.size(); ++i){
			AvalueCopy[i] = Avalue[i] + -2*(*min);
		}
	}
	else{
		for (int i = 0; i < Avalue.size(); ++i){
			AvalueCopy[i] = Avalue[i];
		}
	}
 	for (int i = 0; i < numCol; ++i){
		for (int j = Astart[i]; j < Astart[i + 1]; ++j){
			appendAdj(&adjList[i], Aindex[j] + numCol, AvalueCopy[j]);
			appendAdj(&adjList[Aindex[j] + numCol], i, AvalueCopy[j]); 
		}
	}
} 

// Grab initial partitons (rhs, bounds, obj coeff)
void HModel::initEQs(){
	int numParts = 0;
	initStorage();
	lp2Graph();
	set<double> rhs;
	set<double>::iterator rhsIdx;
	set<tuple<double, double, double > > objBounds;
	set<tuple<double, double, double > >::iterator objBoundsIdx;
	int varColor = 0;
	int conColor = numCol;
	// Prefill vector 
	for (int i = 0; i < numRow + numCol; ++i){
		initialParts.push_back(0);
	}
	// Grab right hand sides of constraints
	for (int i = 0; i < numRow; ++i){
		rhs.insert(Rhs[i]);
	}
	// Tuples of obj coeff, lower and upper bounds
	for (int i = 0; i < numCol; ++i){
		objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i]));
	}
	// Assign initial var colors
	for (objBoundsIdx = objBounds.begin(); objBoundsIdx != objBounds.end(); ++objBoundsIdx){
		for (int i = 0; i < numCol; ++i){
			if (colCost[i] == get<0>(*objBoundsIdx) && colLower[i] == 
				get<1>(*objBoundsIdx) && colUpper[i] == get<2>(*objBoundsIdx)){
				initialParts[i] = varColor;
			}
		}
		S.push(varColor);
		SCheck[varColor] = true;
		varColor++;
	}
	// Assign initial con colors;
	for (rhsIdx = rhs.begin(); rhsIdx != rhs.end(); ++rhsIdx){
		for (int i = 0; i < numRow; ++i){
			if (Rhs[i] == *rhsIdx){
				initialParts[i + numCol] = conColor;
			}
		}
		S.push(conColor);
		SCheck[conColor] = true;
		conColor++;
	}
	// Define number of colors used so far
	vCol = varColor;
	cCol = conColor;
	for (int i = 0; i < numRow + numCol; ++i){
		append(&C[initialParts[i]], i);
		color.push_back(initialParts[i]);
	}
}

// Discretize partitions
void HModel::computeEQs(){
	/* 41 - 44 -- Note: Stack is already sorted because
	of how it is created */ 
	colorsAdj = NULL;
	numParts = 0;
	/* 45 - */
	while (!S.empty()){
		r = S.top();
		S.pop();
		SCheck[r] = false;
		for (v = C[r]; v != NULL; v = v->next){
			for (w = adjList[v->data]; w != NULL; w = w->next){
				cdeg[w->label] += w->w;
				isAdj[w->label]++;
				if (isAdj[w->label] == 1){
					A[color[w->label]]->push_back(w->label);
				}
				if (!exists(&colorsAdj, color[w->label])){
					append(&colorsAdj, color[w->label]);
				}
				if (cdeg[w->label] > maxcdeg[color[w->label]]){
					maxcdeg[color[w->label]] = cdeg[w->label];
				}
			}
		}
		for (c = colorsAdj; c != NULL; c = c->next){
			if (listSize(&C[c->data]) != A[c->data]->size()){
				mincdeg[c->data] = 0;
			}
			else{
				mincdeg[c->data] = maxcdeg[c->data];
				for (u = A[c->data]->begin(); u != A[c->data]->end(); ++u){
					if (cdeg[*u] < mincdeg[c->data]) mincdeg[c->data] = cdeg[*u];
				}
			}
		}
		list<int> colorsSplit;
		for (c = colorsAdj; c!= NULL; c = c->next){
			if (mincdeg[c->data] < maxcdeg[c->data]){
				colorsSplit.push_back(c->data);
			}
		}
		colorsSplit.sort();
		for (s = colorsSplit.begin(); s != colorsSplit.end(); ++s){
			splitColor(*s);
		}
		// Reset attributes
		for (c = colorsAdj; c != NULL; c = c->next){
			for (u = A[c->data]->begin(); u != A[c->data]->end(); ++u){
				cdeg[*u] = 0;
				isAdj[*u] = 0;
			}
			maxcdeg[c->data] = 0;
			A[c->data]->clear();
			deleteNode(&colorsAdj, c);
		}
		for (int i = 0; i < C.size(); ++i)
			if (C[i] != NULL && listSize(&C[i]) == 1) isolates[C[i]->data] = true;
	} 
	// for (int i = 0; i < C.size(); ++i){
	// 	if (C[i] != NULL){
	// 		cout << "color: " << i << endl;
	// 		printList(&C[i]);
	// 		cout << "\n" << endl;
	// 	}
	// }
	aggClear();
	for (int i = 0; i < C.size(); ++i){
		if (C[i] != NULL){
			if (i >= numCol){
				conColorReps.push_back(C[i]->data);
				aggNumRow ++;
				aggNumTot ++;
				aggRowIdx.push_back(i);
			}
			else{
				reps.push_back(C[i]->data);
				aggNumCol ++;
				aggNumTot ++;
				aggColIdx.push_back(i);
				for (v = C[i]; v != NULL; v = v->next){
					vColor[v->data] = i;
				}
			}
		}
	}
	getNewRows();
	aggregate();
}

// Color splitting (cell split) subroutine
void HModel::splitColor(int s){
	set<double> cdegCopy;
	vector<int> colorFreq(numRow + numCol, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	bool varOrCon = (s < numCol) ? true : false;
	degSumColor.insert(pair<double, int>(mincdeg[s], s)); 
	colorFreq[0] = listSize(&C[s]) - A[s]->size(); 
	for (u = A[s]->begin(); u != A[s]->end(); ++u){
		if (varOrCon){
			ret = degSumColor.insert(pair<double, int>(cdeg[*u], vCol));
			ret.second ? vCol++ : vCol;
		}
		else{
			ret = degSumColor.insert(pair<double, int>(cdeg[*u], cCol));
			ret.second ? cCol++ : cCol;
		}
	}
	for (u = A[s]->begin(); u != A[s]->end(); ++u){
			colorFreq[degSumColor[cdeg[*u]]]++;
	}
	int b = distance(colorFreq.begin(), max_element(colorFreq.begin(), colorFreq.end()));
	int instack = (SCheck[s]) ? 1 : 0;
	for (map<double, int>::iterator it = degSumColor.begin(); it != degSumColor.end(); ++it){
		if (it->first == mincdeg[s]){
			if (!instack && it->second != b){
				S.push(it->second);
				SCheck[s] = true;
			}
		}
		else{
			if (instack || it->second != b){
				S.push(it->second);
				SCheck[s] = true;
			}
		}
	}
	for (u = A[s]->begin(); u != A[s]->end(); ++u){
		if (degSumColor[cdeg[*u]] != s){
			del = findNode(&C[s], *u);
			deleteNode(&C[s], del);
			append(&C[degSumColor[cdeg[*u]]], *u);
			color[*u] = degSumColor[cdeg[*u]];
		}
	}
}

// Isolation function
void HModel::isolate(int i){
	del = findNode(&C[color[i]], i);
	int newCol = vCol;
	vCol++;
	deleteNode(&C[color[i]], del);
	append(&C[newCol], i);
	color[i] = newCol;
	for (int i = 0; i < C.size(); ++i){
		if (C[i] != NULL){
			if (listSize(&C[i]) == 1){
				isolates[C[i]->data] = true;
			}
		}
	}
	SCheck[newCol] = true;
	S.push(newCol);
}

// clear out storage for aggregated model
void HModel::aggClear(){
	aggNumRow = 0;
	aggNumCol = 0;
	aggNumTot = 0;
	numNewRows = 0;
	numNewCols = 0;
	targ = -1;
	vColor.assign(numCol, -1);
	reps.clear();
	conColorReps.clear();
	aggAdjList.clear();
	aggRowIdx.clear();
	aggColIdx.clear();
	aggAstart.clear();
	aggAvalue.clear();
	aggAindex.clear();
	aggColCost.clear();
 	aggColLower.clear();
    aggColUpper.clear();
    aggRowLower.clear();
    aggRowUpper.clear();
}

// Find coefficients for aggregated variables
vector<int> HModel::getCoeff(int color){
	vector<int> coeff;
	int count = 0;
	for (int i = 0; i < aggRowIdx.size(); ++i){
		for (w = adjList[C[aggRowIdx[i]]->data]; w != NULL; w = w->next){
			if (vColor[w->label] == color)
				count += w->w;
		}
		coeff.push_back(count);
		count = 0;
	}
	for (int i = 0; i < aggAdjList.size(); ++i){
		for (w = aggAdjList[i]; w != NULL; w = w->next){
			if (w->label == color)
				count += w->w;
		}
		coeff.push_back(count);
		count = 0;
	}
	return coeff;
}

// Find obj coefficients for aggregated variables
int HModel::getObj(int color){
	int obj = 0;
	for (int i = 0; i < colCost.size(); ++i){
		if (vColor[i] == color && colCost[i]){
			obj += colCost[i];
		}
	}
	return obj;
}

// Aggregate variables based on current EP
void HModel::aggregate(){
	aggAstart.push_back(0);
	for (int i = 0; i < aggColIdx.size(); ++i){
		vector<int> coeff = getCoeff(aggColIdx[i]);
		for (int j = 0; j < coeff.size(); ++j){
			if (coeff[j]){
				aggAvalue.push_back(coeff[j]);
				aggAindex.push_back(j);
			}
		}
		if (aggColIdx[i] < numCol){
			if (masterIter && prevBasicColor[oldColor[reps[i]]]){
				startingBasis.push_back(i);
				startingBasicValue.push_back(prevBasicValue[oldColor[reps[i]]]);
			}
			aggColLower.push_back(colLower[C[aggColIdx[i]]->data]);
			aggColUpper.push_back(colUpper[C[aggColIdx[i]]->data]);
		}
		else{
			aggColLower.push_back(0);
			aggColUpper.push_back(0);
		}
		aggAstart.push_back(aggAindex.size());
	}
	for (int i = 0; i < aggRowIdx.size() + numNewRows; ++i){
		if (i < aggRowIdx.size()){
			if (masterIter && prevBasicColor[oldColor[conColorReps[i]]]){
				startingBasis.push_back(i + aggNumCol);
				startingBasicValue.push_back(prevBasicValue[oldColor[conColorReps[i]]]);
			}
			aggRowLower.push_back(rowLower[C[aggRowIdx[i]]->data - numCol]);
			aggRowUpper.push_back(rowUpper[C[aggRowIdx[i]]->data - numCol]);
		}
		else{
			aggRowLower.push_back(0);
			aggRowUpper.push_back(0);
		}
	}
	//if (!masterIter){
	for (int i = 0; i < aggColIdx.size(); ++i){
		if (aggColIdx[i] < numCol)
			aggColCost.push_back(getObj(aggColIdx[i]));
		else
			aggColCost.push_back(0);
	}
	// for (int i = 0; i < startingBasis.size(); ++i)
	// 	cout << "var: " << startingBasis[i] << " has basic value: " << startingBasicValue[i] << endl;

	// cout << "matrix" << endl;
	// for (int i = 0; i < aggColIdx.size(); ++i){
	// 	cout << "col: " << aggColIdx[i] << endl;
	// 	for (int j = aggAstart[i]; j < aggAstart[i + 1]; ++j){
	// 		cout << "coeff: " << aggAvalue[j] << " at row: " << aggAindex[j] << endl;
	// 	}
	// }
}

// Check if parition is discretized
bool HModel::discrete(){
	for (int i = 0; i < isolates.size(); ++i){
		if (!isolates[i]){
			iso = i;
			prevColor = color[i];
			oldColor = color;
			return false;
		}
	}
	return true;
}

// Create pairs for residual variables
void HModel::getNewRows(){
	if (!masterIter)
		return;
	vector<int> leftOvers;
	for (int i = 0; i < reps.size(); ++i){
		if (oldColor[reps[i]] == oldColor[iso] && reps[i] != iso){
			aggAdjList.push_back(NULL);
			appendAdj(&aggAdjList[numNewRows], color[iso], 1);
			appendAdj(&aggAdjList[numNewRows], color[reps[i]], -1);
			appendAdj(&aggAdjList[numNewRows], numTot + numNewCols, -1);
			aggColIdx.push_back(numTot + numNewCols - 1);
			numNewRows++;
			numNewCols++;
		}
		else if (oldColor[reps[i]] != oldColor[iso] && reps[i] != iso){
			leftOvers.push_back(reps[i]);
		}
	}
	sort(leftOvers.begin(), leftOvers.end(), greater<int>());
 	while(!leftOvers.empty()){
		targ = leftOvers.back();
		leftOvers.pop_back();
		for (int i = 0; i < leftOvers.size(); ++i){
			if (oldColor[targ] == oldColor[leftOvers[i]]){
				aggAdjList.push_back(NULL);
				appendAdj(&aggAdjList[numNewRows], color[targ], 1);
				appendAdj(&aggAdjList[numNewRows], color[leftOvers[i]], -1);
				appendAdj(&aggAdjList[numNewRows], numTot + numNewCols, -1);
				aggColIdx.push_back(numTot + numNewCols - 1);
				numNewRows++;
				numNewCols++;
			}
		}
	}
	aggNumTot += numNewRows + numNewCols;
	aggNumCol += numNewCols;
	aggNumRow += numNewRows;
}

// EP until discrete partition
void HModel::equitable(){
	computeEQs();
	// isoIter = 0;
	// for (int i = 0; i < isolates.size(); ++i){
	// 	if (isolates[i] == false){
	// 		isoIter ++;
	// 		isolate(i);
	// 		cin.get();
	// 	}
	// }
}

void HModel::setup_transposeLP() {
    if (intOption[INTOPT_TRANSPOSE_FLAG] == 0)
        return;

    int transposeCancelled = 0;
    if (1.0 * aggNumCol / aggNumRow > 0.2) {
//        cout << "transpose-cancelled-by-ratio" << endl;
        transposeCancelled = 1;
        return;
    }

    // Convert primal cost to dual bound
    const double inf = HSOL_CONST_INF;
    vector<double> dualRowLower(aggNumCol);
    vector<double> dualRowUpper(aggNumCol);
    for (int j = 0; j < aggNumCol; j++) {
        double lower = aggColLower[j];
        double upper = aggColUpper[j];

        /*
         * Primal      Dual
         * Free        row = c
         * x > 0       row < c
         * x < 0       row > c
         * x = 0       row free
         * other       cancel
         */

        if (lower == -inf && upper == inf) {
            dualRowLower[j] = aggColCost[j];
            dualRowUpper[j] = aggColCost[j];
        } else if (lower == 0 && upper == inf) {
            dualRowLower[j] = -inf;
            dualRowUpper[j] = aggColCost[j];
        } else if (lower == -inf && upper == 0) {
            dualRowLower[j] = aggColCost[j];
            dualRowUpper[j] = +inf;
        } else if (lower == 0 && upper == 0) {
            dualRowLower[j] = -inf;
            dualRowUpper[j] = +inf;
        } else {
            transposeCancelled = 1;
            break;
        }
    }

    // Check flag
    if (transposeCancelled == 1) {
//        cout << "transpose-cancelled-by-column" << endl;
        return;
    }

    // Convert primal row bound to dual variable cost
    vector<double> dualColLower(aggNumCol);
    vector<double> dualColUpper(aggNumCol);
    vector<double> dualCost(aggNumCol);
    for (int i = 0; i < aggNumCol; i++) {
        double lower = aggRowLower[i];
        double upper = aggRowUpper[i];

        /*
         * Primal      Dual
         * row = b     Free
         * row < b     y < 0
         * row > b     y > 0
         * row free    y = 0
         * other       cancel
         */

        if (lower == upper) {
            dualColLower[i] = -inf;
            dualColUpper[i] = +inf;
            dualCost[i] = -lower;
        } else if (lower == -inf && upper != inf) {
            dualColLower[i] = -inf;
            dualColUpper[i] = 0;
            dualCost[i] = -upper;
        } else if (lower != -inf && upper == inf) {
            dualColLower[i] = 0;
            dualColUpper[i] = +inf;
            dualCost[i] = -lower;
        } else if (lower == -inf && upper == inf) {
            dualColLower[i] = 0;
            dualColUpper[i] = 0;
            dualCost[i] = 0;
        } else {
            transposeCancelled = 1;
            break;
        }
    }

    // Check flag
    if (transposeCancelled == 1) {
//        cout << "transpose-cancelled-by-row" << endl;
        return;
    }

    // We can now really transpose things
    vector<int> iwork(aggNumCol, 0);
    vector<int> ARstart(aggNumCol + 1, 0);
    int AcountX = aggAindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++)
        iwork[aggAindex[k]]++;
    for (int i = 1; i <= aggNumCol; i++)
        ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < aggNumCol; i++)
        iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < aggNumCol; iCol++) {
        for (int k = aggAstart[iCol]; k < aggAstart[iCol + 1]; k++) {
            int iRow = aggAindex[k];
            int iPut = iwork[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = aggAvalue[k];
        }
    }

    // Transpose the problem!
    swap(aggNumCol, aggNumCol);
    aggAstart.swap(ARstart);
    aggAindex.swap(ARindex);
    aggAvalue.swap(ARvalue);
    aggColLower.swap(dualColLower);
    aggColUpper.swap(dualColUpper);
    aggRowLower.swap(dualRowLower);
    aggRowUpper.swap(dualRowUpper);
    aggColCost.swap(dualCost);
//    cout << "problem-transposed" << endl;
}

void HModel::setup_scaleMatrix() {
    if (intOption[INTOPT_SCALE_FLAG] == 0)
        return;

    // Reset all scaling to 1
    vector<double> colScale(aggNumCol, 1);
    vector<double> rowScale(aggNumRow, 1);

    // Find out min0 / max0, skip on if in [0.2, 5]
    const double inf = HSOL_CONST_INF;
    double min0 = inf, max0 = 0;
    for (int k = 0, AnX = aggAstart[aggNumCol]; k < AnX; k++) {
        double value = fabs(aggAvalue[k]);
        min0 = min(min0, value);
        max0 = max(max0, value);
    }
    if (min0 >= 0.2 && max0 <= 5)
        return;

    // See if we want to include cost include if min-cost < 0.1
    double minc = inf;
    for (int i = 0; i < aggNumCol; i++)
        if (aggColCost[i])
            minc = min(minc, fabs(aggColCost[i]));
    bool doCost = minc < 0.1;

    // Search up to 6 times
    vector<double> rowMin(aggNumRow, inf);
    vector<double> rowMax(aggNumRow, 1 / inf);
    for (int search_count = 0; search_count < 6; search_count++) {
        // Find column scale, prepare row data
        for (int iCol = 0; iCol < aggNumCol; iCol++) {
            // For column scale (find)
            double colMin = inf;
            double colMax = 1 / inf;
            double myCost = fabs(aggColCost[iCol]);
            if (doCost && myCost != 0)
                colMin = min(colMin, myCost), colMax = max(colMax, myCost);
            for (int k = aggAstart[iCol]; k < aggAstart[iCol + 1]; k++) {
                double value = fabs(aggAvalue[k]) * rowScale[aggAindex[k]];
                colMin = min(colMin, value), colMax = max(colMax, value);
            }
            colScale[iCol] = 1 / sqrt(colMin * colMax);

            // For row scale (only collect)
            for (int k = aggAstart[iCol]; k < aggAstart[iCol + 1]; k++) {
                int iRow = aggAindex[k];
                double value = fabs(aggAvalue[k]) * colScale[iCol];
                rowMin[iRow] = min(rowMin[iRow], value);
                rowMax[iRow] = max(rowMax[iRow], value);
            }
        }

        // For row scale (find)
        for (int iRow = 0; iRow < aggNumRow; iRow++)
            rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
        rowMin.assign(aggNumRow, inf);
        rowMax.assign(aggNumRow, 1 / inf);
    }

    // Make it numerical better
    const double ln2 = log(2.0);
    for (int iCol = 0; iCol < aggNumCol; iCol++)
        colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
    for (int iRow = 0; iRow < aggNumRow; iRow++)
        rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));

    // Apply scaling to matrix and bounds
    for (int iCol = 0; iCol < aggNumCol; iCol++)
        for (int k = aggAstart[iCol]; k < aggAstart[iCol + 1]; k++)
            aggAvalue[k] *= (colScale[iCol] * rowScale[aggAindex[k]]);

    for (int iCol = 0; iCol < aggNumCol; iCol++) {
        aggColLower[iCol] /= aggColLower[iCol] == -inf ? 1 : colScale[iCol];
        aggColUpper[iCol] /= aggColUpper[iCol] == +inf ? 1 : colScale[iCol];
        aggColCost[iCol] *= colScale[iCol];
    }
    for (int iRow = 0; iRow < aggNumRow; iRow++) {
        aggRowLower[iRow] *= aggRowLower[iRow] == -inf ? 1 : rowScale[iRow];
        aggRowUpper[iRow] *= aggRowUpper[iRow] == +inf ? 1 : rowScale[iRow];
    }
}

void HModel::setup_tightenBound() {
    if (intOption[INTOPT_TIGHT_FLAG] == 0)
        return;

    // Make a AR copy
    vector<int> iwork(aggNumCol, 0);
    vector<int> ARstart(aggNumCol + 1, 0);
    int AcountX = aggAindex.size();
    vector<int> ARindex(AcountX);
    vector<double> ARvalue(AcountX);
    for (int k = 0; k < AcountX; k++)
        iwork[aggAindex[k]]++;
    for (int i = 1; i <= aggNumCol; i++)
        ARstart[i] = ARstart[i - 1] + iwork[i - 1];
    for (int i = 0; i < aggNumCol; i++)
        iwork[i] = ARstart[i];
    for (int iCol = 0; iCol < aggNumCol; iCol++) {
        for (int k = aggAstart[iCol]; k < aggAstart[iCol + 1]; k++) {
            int iRow = aggAindex[k];
            int iPut = iwork[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = aggAvalue[k];
        }
    }

    // Save column bounds
    vector<double> colLower0 = aggColLower;
    vector<double> colUpper0 = aggColUpper;

    double big_B = 1e10;
    int iPass = 0;
    for (;;) {
        int numberChanged = 0;
        for (int iRow = 0; iRow < aggNumCol; iRow++) {
            // SKIP free rows
            if (aggRowLower[iRow] < -big_B && aggRowUpper[iRow] > big_B)
                continue;

            // possible row
            int ninfU = 0;
            int ninfL = 0;
            double xmaxU = 0.0;
            double xminL = 0.0;
            int myStart = ARstart[iRow];
            int myEnd = ARstart[iRow + 1];
            // Compute possible lower and upper ranges

            for (int k = myStart; k < myEnd; ++k) {
                int iCol = ARindex[k];
                double value = ARvalue[k];
                double upper = value > 0 ? aggColUpper[iCol] : -aggColLower[iCol];
                double lower = value > 0 ? aggColLower[iCol] : -aggColUpper[iCol];
                value = fabs(value);
                if (upper < big_B)
                    xmaxU += upper * value;
                else
                    ++ninfU;
                if (lower > -big_B)
                    xminL += lower * value;
                else
                    ++ninfL;
            }

            // Build in a margin of error
            xmaxU += 1.0e-8 * fabs(xmaxU);
            xminL -= 1.0e-8 * fabs(xminL);

            double xminLmargin =
                    (fabs(xminL) > 1.0e8) ? 1e-12 * fabs(xminL) : 0;
            double xmaxUmargin =
                    (fabs(xmaxU) > 1.0e8) ? 1e-12 * fabs(xmaxU) : 0;

            // Skip redundant row : also need to consider U < L  case
            double comp_U = xmaxU + ninfU * 1.0e31;
            double comp_L = xminL - ninfL * 1.0e31;
            if (comp_U <= aggRowUpper[iRow] + 1e-7
                    && comp_L >= aggRowLower[iRow] - 1e-7)
                continue;

            double row_L = aggRowLower[iRow];
            double row_U = aggRowUpper[iRow];

            // Now see if we can tighten column bounds
            for (int k = myStart; k < myEnd; ++k) {
                double value = ARvalue[k];
                int iCol = ARindex[k];
                double col_L = aggColLower[iCol];
                double col_U = aggColUpper[iCol];
                double new_L = -HSOL_CONST_INF;
                double new_U = +HSOL_CONST_INF;

                if (value > 0.0) {
                    if (row_L > -big_B && ninfU <= 1
                            && (ninfU == 0 || col_U > +big_B))
                        new_L = (row_L - xmaxU) / value + (1 - ninfU) * col_U
                                - xmaxUmargin;
                    if (row_U < +big_B && ninfL <= 1
                            && (ninfL == 0 || col_L < -big_B))
                        new_U = (row_U - xminL) / value + (1 - ninfL) * col_L
                                + xminLmargin;
                } else {
                    if (row_L > -big_B && ninfU <= 1
                            && (ninfU == 0 || col_L < -big_B))
                        new_U = (row_L - xmaxU) / value + (1 - ninfU) * col_L
                                + xmaxUmargin;
                    if (row_U < +big_B && ninfL <= 1
                            && (ninfL == 0 || col_U > +big_B))
                        new_L = (row_U - xminL) / value + (1 - ninfL) * col_U
                                - xminLmargin;
                }

                if (new_U < col_U - 1.0e-12 && new_U < big_B) {
                    aggColUpper[iCol] = max(new_U, col_L);
                    numberChanged++;
                }
                if (new_L > col_L + 1.0e-12 && new_L > -big_B) {
                    aggColLower[iCol] = min(new_L, col_U);
                    numberChanged++;
                }
            }
        }

        if (numberChanged == 0)
            break;
        iPass++;
        if (iPass > 10)
            break;

    }

    double useTolerance = 1.0e-3;
    for (int iCol = 0; iCol < aggNumCol; iCol++) {
        if (colUpper0[iCol] > colLower0[iCol] + useTolerance) {
            const double relax = 100.0 * useTolerance;
            if (aggColUpper[iCol] - aggColLower[iCol] < useTolerance + 1.0e-8) {
                aggColLower[iCol] = max(colLower0[iCol], aggColLower[iCol] - relax);
                aggColUpper[iCol] = min(colUpper0[iCol], aggColUpper[iCol] + relax);
            } else {
                if (aggColUpper[iCol] < colUpper0[iCol]) {
                    aggColUpper[iCol] = min(aggColUpper[iCol] + relax,
                            colUpper0[iCol]);
                }
                if (aggColLower[iCol] > colLower0[iCol]) {
                    aggColLower[iCol] = min(aggColLower[iCol] - relax,
                            colLower0[iCol]);
                }
            }
        }
    }
}

void HModel::setup_shuffleColumn() {
    if (intOption[INTOPT_PERMUTE_FLAG] == 0)
        return;

    // 1. Shuffle the column index
    HRandom localRandom;
    for (int i = 0; i < 10; i++)
        localRandom.intRandom();
    vector<int> iFrom(aggNumCol);
    for (int i = 0; i < aggNumCol; i++)
        iFrom[i] = i;
    for (int i = aggNumCol - 1; i >= 1; i--) {
        int j = localRandom.intRandom() % (i + 1);
        swap(iFrom[i], iFrom[j]);
    }

    // 2. Save original copy
    vector<int> start = aggAstart;
    vector<int> index = aggAindex;
    vector<double> value = aggAvalue;
    vector<double> lower = aggColLower;
    vector<double> upper = aggColUpper;
    vector<double> xcost = aggColCost;
    vector<int> ibreak = intBreak;
    vector<double> dxpert = dblXpert;

    // 3. Generate the permuted matrix
    int countX = 0;
    for (int i = 0; i < aggNumCol; i++) {
        int ifrom = iFrom[i];
        aggAstart[i] = countX;
        for (int k = start[ifrom]; k < start[ifrom + 1]; k++) {
            aggAindex[countX] = index[k];
            aggAvalue[countX] = value[k];
            countX++;
        }
        aggColLower[i] = lower[ifrom];
        aggColUpper[i] = upper[ifrom];
        aggColCost[i] = xcost[ifrom];
        intBreak[i] = ibreak[ifrom];
        dblXpert[i] = dxpert[ifrom];
    }
    assert(aggAstart[aggNumCol] == countX);
}

void HModel::setup_allocWorking() {
    // Setup starting base
    basicIndex.resize(aggNumRow);
    for (int iRow = 0; iRow < aggNumRow; iRow++)
        basicIndex[iRow] = iRow + aggNumCol;
    nonbasicFlag.assign(aggNumTot, 0);
    nonbasicMove.resize(aggNumTot);
    for (int i = 0; i < aggNumCol; i++)
        nonbasicFlag[i] = 1;

    // Matrix, factor
    matrix.setup(aggNumCol, aggNumRow, &aggAstart[0], &aggAindex[0], &aggAvalue[0]);
    factor.setup(aggNumCol, aggNumRow, &aggAstart[0], &aggAindex[0], &aggAvalue[0],
            &basicIndex[0]);
    limitUpdate = 5000;

    // Setup other buffer
    buffer.setup(aggNumRow);
    bufferLong.setup(aggNumCol);

    // Setup bounds and solution spaces
    workCost.assign(aggNumTot, 0);
    workDual.assign(aggNumTot, 0);
    workShift.assign(aggNumTot, 0);

    workLower.assign(aggNumTot, 0);
    workUpper.assign(aggNumTot, 0);
    workRange.assign(aggNumTot, 0);
    workValue.assign(aggNumTot, 0);

    baseLower.assign(aggNumRow, 0);
    baseUpper.assign(aggNumRow, 0);
    baseValue.assign(aggNumRow, 0);
}

void HModel::initCost(int perturb) {
    // Copy the cost
    for (int i = 0; i < aggNumCol; i++)
        workCost[i] = aggColCost[i];
    for (int i = aggNumCol; i < aggNumTot; i++)
        workCost[i] = 0;
    workShift.assign(aggNumTot, 0);

    // See if we want to skip perturbation
    problemPerturbed = 0;
    if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0)
        return;
    problemPerturbed = 1;

    // Perturb the original costs, scale down if is too big
    double bigc = 0;
    for (int i = 0; i < aggNumCol; i++)
        bigc = max(bigc, fabs(workCost[i]));
    if (bigc > 100)
        bigc = sqrt(sqrt(bigc));

    // If there's few boxed variables, we will just use Simple perturbation
    double boxedRate = 0;
    for (int i = 0; i < aggNumTot; i++)
        boxedRate += (workRange[i] < 1e30);
    boxedRate /= aggNumTot;
    if (boxedRate < 0.01)
        bigc = min(bigc, 1.0);
    if (bigc < 1) {
//        bigc = sqrt(bigc);
    }

    // Determine the perturbation base
    double base = 5e-7 * bigc;

    // Now do the perturbation
    for (int i = 0; i < aggNumCol; i++) {
        double lower = aggColLower[i];
        double upper = aggColUpper[i];
        double xpert = (fabs(workCost[i]) + 1) * base * (1 + dblXpert[i]);
        if (lower == -HSOL_CONST_INF && upper == HSOL_CONST_INF) {
            // Free - no perturb
        } else if (upper == HSOL_CONST_INF) {  // Lower
            workCost[i] += xpert;
        } else if (lower == -HSOL_CONST_INF) { // Upper
            workCost[i] += -xpert;
        } else if (lower != upper) {  // Boxed
            workCost[i] += (workCost[i] >= 0) ? xpert : -xpert;
        } else {
            // Fixed - no perturb
        }
    }

    for (int i = aggNumCol; i < aggNumTot; i++) {
        workCost[i] += (0.5 - dblXpert[i]) * 1e-12;
    }
}

void HModel::initBound(int phase) {
    // Copy bounds
    for (int i = 0; i < aggNumCol; i++) {
        workLower[i] = aggColLower[i];
        workUpper[i] = aggColUpper[i];
    }
    for (int i = 0, j = aggNumCol; i < aggNumRow; i++, j++) {
        workLower[j] = -aggRowUpper[i];
        workUpper[j] = -aggRowLower[i];
    }

    // Change to dual phase 1 bound
    if (phase == 1) {
        const double inf = HSOL_CONST_INF;
        for (int i = 0; i < aggNumTot; i++) {
            if (workLower[i] == -inf && workUpper[i] == inf) {
                // Won't change for row variables: they should never
                // Become non basic
                if (i >= aggNumCol)
                    continue;
                workLower[i] = -1000, workUpper[i] = 1000;  // FREE
            } else if (workLower[i] == -inf) {
                workLower[i] = -1, workUpper[i] = 0;        // UPPER
            } else if (workUpper[i] == inf) {
                workLower[i] = 0, workUpper[i] = 1;         // LOWER
            } else {
                workLower[i] = 0, workUpper[i] = 0;         // BOXED or FIXED
            }
        }
    }

    // Work out ranges
    for (int i = 0; i < aggNumTot; i++)
        workRange[i] = workUpper[i] - workLower[i];
}

void HModel::initValue() {
    for (int i = 0; i < aggNumTot; i++) {
        if (nonbasicFlag[i]) {
            if (workLower[i] == workUpper[i]) {
                // Fixed
                workValue[i] = workLower[i];
                nonbasicMove[i] = 0;
            } else if (workLower[i] != -HSOL_CONST_INF) {
                // Lower
                workValue[i] = workLower[i];
                nonbasicMove[i] = 1;
            } else if (workUpper[i] != +HSOL_CONST_INF) {
                // UPPER
                workValue[i] = workUpper[i];
                nonbasicMove[i] = -1;
            } else {
                // FREE
                workValue[i] = 0;
                nonbasicMove[i] = 0;
            }
        } else {
            nonbasicMove[i] = 0;
        }
    }
}

void HModel::computeFactor() {
    try {
        factor.build();
    } catch (runtime_error& error) {
        cout << error.what() << endl;
        problemStatus = 3;
        // TODO Change to return
        writePivots("failed");
        exit(0);
    }
    countUpdate = 0;
}

void HModel::computeDual() {
    buffer.clear();
    for (int iRow = 0; iRow < aggNumRow; iRow++) {
        buffer.index[iRow] = iRow;
        buffer.array[iRow] = workCost[basicIndex[iRow]]
                + workShift[basicIndex[iRow]];
    }
    buffer.count = aggNumRow;
    factor.btran(buffer, 1);

    bufferLong.clear();
    matrix.price_by_col(bufferLong, buffer);
    for (int i = 0; i < aggNumCol; i++)
        workDual[i] = workCost[i] - bufferLong.array[i];
    for (int i = aggNumCol; i < aggNumTot; i++)
        workDual[i] = workCost[i] - buffer.array[i - aggNumCol];
}

void HModel::computeDualInfeasInDual(int *dualInfeasCount) {
    int workCount = 0;
    const double inf = HSOL_CONST_INF;
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    for (int i = 0; i < aggNumTot; i++) {
        // Only for non basic variables
        if (!nonbasicFlag[i])
            continue;
        // Free
        if (workLower[i] == -inf && workUpper[i] == inf)
            workCount += (fabs(workDual[i]) >= tau_d);
        // In dual, assuming that boxed variables will be flipped
        if (workLower[i] == -inf || workUpper[i] == inf)
            workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
    }
    *dualInfeasCount = workCount;
}

void HModel::computeDualInfeasInPrimal(int *dualInfeasCount) {
    int workCount = 0;
    const double inf = HSOL_CONST_INF;
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    for (int i = 0; i < aggNumTot; i++) {
        // Only for non basic variables
        if (!nonbasicFlag[i])
            continue;
        // Free
        if (workLower[i] == -inf && workUpper[i] == inf)
            workCount += (fabs(workDual[i]) >= tau_d);
        // In primal don't assume flip
        workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
    }
    *dualInfeasCount = workCount;
}

void HModel::correctDual(int *freeInfeasCount) {
    const double tau_d = dblOption[DBLOPT_DUAL_TOL];
    const double inf = HSOL_CONST_INF;
    int workCount = 0;
    for (int i = 0; i < aggNumTot; i++) {
        if (nonbasicFlag[i]) {
            if (workLower[i] == -inf && workUpper[i] == inf) {
                // FREE variable
                workCount += (fabs(workDual[i]) >= tau_d);
            } else if (nonbasicMove[i] * workDual[i] <= -tau_d) {
                if (workLower[i] != -inf && workUpper[i] != inf) {
                    // Boxed variable = flip
                    flipBound(i);
                } else {
                    // Other variable = shift
                    problemPerturbed = 1;
                    if (nonbasicMove[i] == 1) {
                        double dual = (1 + random.dblRandom()) * tau_d;
                        double shift = dual - workDual[i];
                        workDual[i] = dual;
                        workCost[i] = workCost[i] + shift;
                    } else {
                        double dual = -(1 + random.dblRandom()) * tau_d;
                        double shift = dual - workDual[i];
                        workDual[i] = dual;
                        workCost[i] = workCost[i] + shift;
                    }
                }
            }
        }
    }

    *freeInfeasCount = workCount;
}

void HModel::computePrimal() {
    buffer.clear();
    for (int i = 0; i < aggNumTot; i++)
        if (nonbasicFlag[i] && workValue[i] != 0)
            matrix.collect_aj(buffer, i, workValue[i]);
    factor.ftran(buffer, 1);
    for (int i = 0; i < aggNumRow; i++) {
        int iCol = basicIndex[i];
        baseValue[i] = -buffer.array[i];
        baseLower[i] = workLower[iCol];
        baseUpper[i] = workUpper[iCol];
    }
}

void HModel::computeObject(int phase) {
    objective = 0;
    for (int i = 0; i < aggNumTot; i++)
        if (nonbasicFlag[i])
            objective += workValue[i] * workDual[i];
    if (phase != 1)
        objective -= objOffset;
}

void HModel::shiftCost(int iCol, double amount) {
    problemPerturbed = 1;
    assert(workShift[iCol] == 0);
    workShift[iCol] = amount;
}

void HModel::shiftBack(int iCol) {
    workDual[iCol] -= workShift[iCol];
    workShift[iCol] = 0;
}

void HModel::flipBound(int iCol) {
    const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
    workValue[iCol] = move == 1 ? workLower[iCol] : workUpper[iCol];
}

void HModel::updateFactor(HVector *column, HVector *row_ep, int *iRow,
        int *hint) {
    timer.recordStart(HTICK_UPDATE_FACTOR);
    factor.update(column, row_ep, iRow, hint);
    if (countUpdate >= limitUpdate)
        *hint = 1;
    timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updateMatrix(int columnIn, int columnOut) {
    timer.recordStart(HTICK_UPDATE_FACTOR);
    matrix.update(columnIn, columnOut);
    timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updatePivots(int columnIn, int rowOut, int sourceOut) {
    int columnOut = basicIndex[rowOut];

    // Incoming variable
    basicIndex[rowOut] = columnIn;
    nonbasicFlag[columnIn] = 0;
    nonbasicMove[columnIn] = 0;
    baseLower[rowOut] = workLower[columnIn];
    baseUpper[rowOut] = workUpper[columnIn];

    // Outgoing variable
    nonbasicFlag[columnOut] = 1;
    if (workLower[columnOut] == workUpper[columnOut]) {
        workValue[columnOut] = workLower[columnOut];
        nonbasicMove[columnOut] = 0;
    } else if (sourceOut == -1) {
        workValue[columnOut] = workLower[columnOut];
        nonbasicMove[columnOut] = 1;
    } else {
        workValue[columnOut] = workUpper[columnOut];
        nonbasicMove[columnOut] = -1;
    }

    countUpdate++;
}

void HModel::changeUpdate(int updateMethod) {
    factor.change(updateMethod);
}

void HModel::reportPivots(int columnIn, int columnOut, double alpha) {
    if (columnIn >= 0)
        numberIteration++;
    historyColumnIn.push_back(columnIn);
    historyColumnOut.push_back(columnOut);
    historyAlpha.push_back(alpha);
}

void HModel::reportStatus(int status) {
    problemStatus = status;
}

void HModel::printMessage(const char *message) {
    if (intOption[INTOPT_PRINT_FLAG])
        printf("%s\n", message);
}

void HModel::printObject() {
    if (intOption[INTOPT_PRINT_FLAG])
        printf("%10d  %20.10e\n", numberIteration, objective);
}

void HModel::printResult() {
    if (problemStatus == 0) {
        printf("HsolRpIt - OPTIMAL %16s %20.10e %10d %10.3f\n", modelName.c_str(),
                objective, numberIteration, totalTime);
        printf("grep_HsolRpIt - OPTIMAL,%16s,%20.10e,%10d,%10.3f,0\n", modelName.c_str(),
                objective, numberIteration, totalTime);
    } else {
        printf("HsolRpIt - NOT-OPT status = %d\n", problemStatus);
        printf("grep_HsolRpIt - NOT-OPT,%16s,%20.10e,%10d,%10.3f,%d\n", modelName.c_str(), objective, numberIteration, totalTime, problemStatus);
    }
}

void HModel::printProgress() {
//    return;
//    static double nextReport = 0;
//    double currentTime = timer.getTime();
//    if (currentTime >= nextReport) {
//        computeObject();
//        printf("PROGRESS %16s %20.10e %10d %10.3f\n", modelName.c_str(),
//                objective, numberIteration, timer.getTime());
//        if (currentTime < 50) {
//            nextReport = ((int) (5 * currentTime + 1)) / 5.0 - 0.00001;
//        } else if (currentTime < 500) {
//            nextReport = ((int) (currentTime + 1)) - 0.00001;
//        } else {
//            nextReport = ((int) (0.2 * currentTime + 1)) / 0.2 - 0.00001;
//        }
//    }
}

void HModel::writePivots(const char *suffix) {
    string filename = "z-" + modelName + "-" + suffix;
    ofstream output(filename.c_str());
    int count = historyColumnIn.size();
    output << modelName << " " << count << "\t" << totalTime << endl;
    output << setprecision(12);
    for (int i = 0; i < count; i++) {
        output << historyColumnIn[i] << "\t";
        output << historyColumnOut[i] << "\t";
        output << historyAlpha[i] << endl;
    }
    output.close();
}

void HModel::writeMPS(const char *filename) {

}
