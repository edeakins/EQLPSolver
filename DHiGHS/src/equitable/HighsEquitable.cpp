#include "equitable/HighsEquitable.h"
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




HighsEquitable::HighsEquitable(const HighsLp& lp){
	// Original Lp info but edited for cuts
    nCols = lp.numCol_;
	nRows = lp.numRow_;
    nTot = lp.numCol_ + lp.numRow_;
	//nTotal = nTot;
	nnz = lp.nnz_;
    colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
    colLower.assign(lp.colLower_.begin(), lp.colLower_.end());
    colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
    rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
    rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
    Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
    Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
    Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
	colNames.assign(lp.col_names_.begin(), lp.col_names_.end());
	rowNames.assign(lp.row_names_.begin(), lp.row_names_.end());
	// Make color scheme edits for new cut rows
	// These rows are singletons and therefore do not effect anything 
	// other than that theyu exist.
	AindexP.resize(Aindex.size());
    AvaluePos.resize(Avalue.size());
    ARvaluePos.resize(Avalue.size());
    SCheck.assign(nTot, false);
	mincdeg.assign(nTot, 0);
	maxcdeg.assign(nTot, 0);
	cdeg.assign(nTot, 0);
    isAdj.assign(nTot, false);
    numEdges.assign(nTot, 0);
	isolates.assign(nTot, false);
    color.assign(nTot, 0);
	Csize.assign(nTot, 0);
	Asize.assign(nTot, 0);
    C.resize(nTot);
    A.resize(nTot);
	transpose();
	lp2Graph();
	doSaucyEquitable();
	// colorAlloc();
	handleNegatives();
    initRefinement();
}

void HighsEquitable::lp2Graph(){
	int i, j, nColor = 0;
	double sum;
	// Use to color vertices by rhs, bounds, and obj coeff
	set<tuple<double, double> > rhs;
	set<tuple<double, double, double> > objBounds;
	g = (struct amorph_graph *)malloc(sizeof(struct amorph_graph));
	g->sg.n = nTot;
	g->sg.e = nnz;
	int* adj = (int *)calloc(nTot + 1, sizeof(int));
	int* edg = (int *)calloc(2 * nnz, sizeof(int));
	double* w8t = (double *)calloc(2 * nnz, sizeof(double));
	int* colors = (int *)calloc(nTot, sizeof(int));
	char** var_names = (char** )malloc(nTot*sizeof(char*));
	//marks = (char*)calloc(nTot, sizeof(char));
	g->sg.adj = adj;
	g->sg.edg = edg;
	g->sg.w8t = w8t;
	g->colors = colors;
	g->var_names = var_names;
	g->consumer = amorph_print_automorphism;
	// Fill adj, should be the same as Astart and ARstart + nnz
	for (i = 0; i < nCols; ++i)
		adj[i + 1] = Astart[i + 1];
	for (i = 0; i < nRows; ++i)
		adj[i + nCols + 1] = ARstart[i + 1] + nnz;
	// Fill edg and w8t
	for (i = 0; i < nCols; ++i){
		for (j = Astart[i]; j < Astart[i + 1]; ++j){
			edg[j] = Aindex[j] + nCols;
			w8t[j] = Avalue[j];
		}
	}
	for (i = 0; i < nRows; ++i){
		for (j = ARstart[i]; j < ARstart[i + 1]; ++j){
			edg[j + nnz] = ARindex[j];
			w8t[j + nnz] = ARvalue[j];
		}
	}
	// Fill initial colors
	objBounds.insert(make_tuple(colLower[0], colUpper[0], colCost[0]));
	colors[0] = nColor;
	for (i = 1; i < nCols; ++i){
		if (objBounds.insert(make_tuple(colLower[i], colUpper[i], colCost[i])).second)
			colors[i] = ++nColor;
		else
			colors[i] = nColor;
	}
	rhs.insert(make_tuple(rowLower[0], rowUpper[0]));
	colors[nCols] = ++nColor;
	for (i = 1; i < nRows; ++i){
		if (rhs.insert(make_tuple(rowLower[i], rowUpper[i])).second)
			colors[i + nCols] = ++nColor;
		else
			colors[i + nCols] = nColor;
	}	
	// Fill in var_names
	for (i = 0; i < nCols; ++i){
		int tempk = colNames[i].length();
		var_names[i] = (char*)malloc(tempk*sizeof(char));
		strcpy(var_names[i], colNames[i].c_str());
	}
	for(i = 0; i < nRows; ++i)
    {
        int tempk = rowNames[i].length();
        var_names[i + nCols] = (char *)malloc( tempk*sizeof(char) );
        strcpy(var_names[i + nCols], rowNames[i].c_str());
    }
}

void HighsEquitable::doSaucyEquitable(){
	s = saucy_alloc(nTot);
	saucy_search(s, &g->sg, 0, g->colors, on_automorphism, g, &sstats);
}

void HighsEquitable::initRefinement(){
    int i,j;
	int max = 0;
	int numParts = 0;
	int varColor = 0;
	int conColor = nCols;
	set<double> Rhs;
	vector<double> Rhs_;
	set<tuple<double, double, double> > rhs;
	set<tuple<double, double, double, double> > objBounds;
    initialParts.resize(nTot);

    // Partition the rows
	pair<set<tuple<double, double, double> >::iterator, bool> retRow;
	int rowSum = 0;
	for (i = ARstart[0]; i < ARstart[1]; ++i)
		rowSum += ARvaluePos[i];
	rhs.insert(make_tuple(rowLower[0], rowUpper[0], rowSum));
	initialParts[nCols] = conColor;
	S.push(conColor);
	SCheck[conColor] = true;
	rowSum = 0;
    for (i = 1; i < nRows; ++i){
		for (j = ARstart[i]; j < ARstart[i + 1]; ++j)
			rowSum += ARvaluePos[j];
		retRow = rhs.insert(make_tuple(rowLower[i], rowUpper[i], rowSum));
		if (retRow.second){
			initialParts[i + nCols] = ++conColor;
			S.push(conColor);
			SCheck[conColor] = true;
		}
		else{
			initialParts[i + nCols] = conColor;
		}
		rowSum = 0;
	}
	conColor++;

    // Partition the columns
	pair<set<tuple<double, double, double, double> >::iterator, bool> retCol;
	int colSum = 0;
	for (i = Astart[0]; i < Astart[1]; ++i){
		colSum += AvaluePos[i];
	}
	objBounds.insert(make_tuple(colCost[0], colLower[0], colUpper[0], colSum));
	initialParts[0] = varColor;
	S.push(varColor);
	SCheck[varColor] = true;
	colSum = 0;
    for (i = 1; i < nCols; ++i){
		for (j = Astart[i]; j < Astart[i + 1]; ++j)
			colSum += AvaluePos[j];
		retCol = objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i], colSum));
		if (retCol.second){
			initialParts[i] = ++varColor;
			S.push(varColor);
			SCheck[varColor] = true;;
		}
		else{
			initialParts[i] = varColor;
		}
		colSum = 0;
	}
	varColor++;

    // Define number of colors that have been used so far
	vCol = varColor;
	cCol = conColor;

	// Store initial refinement
	for (i = 0; i < nTot; ++i){
		C[initialParts[i]].push_back(i);
		Csize[initialParts[i]]++;
		color[i] = initialParts[i];
	}
	// /*  Count the initial frequency of each color class */	
	// for (i = 0; i < nTot; ++i){
	// 	ccount[color[i]]++;
	// 	if (max < color[i]) max = color[i];
	// }
	// /* Build cell lengths */
	// coloring.clen[0] = ccount[0] - 1;
	// for (i = 0; i < max; ++i) {
	// 	coloring.clen[ccount[i]] = ccount[i+1] - 1;
	// 	ccount[i+1] += ccount[i];
	// }
	// /* Build the label array */
	// for (i = 0; i < nTot; ++i) {
	// 	setLabel(--ccount[color[i]], i);
	// }
	// /* Clear out ccount */
	// for (i = 0; i <= max; ++i) {
	// 	ccount[i] = 0;
	// }
	// /* Update refinement stuff based on initial partition */
	// for (i = 0; i < nTot; i += coloring.clen[i]+1) {
	// 	addInduce(i);
	// 	fixFronts(i, i);
	// }
	// /* Prepare lists based on cell lengths */
	// for (i = 0, j = -1; i < nTot; i += coloring.clen[i] + 1) {
	// 	if (!coloring.clen[i]) continue;
	// 	prevnon[i] = j;
	// 	nextnon[j] = i;
	// 	j = i;
	// }
	// /* Fix the end */
	// prevnon[nTot] = j;
	// nextnon[j] = nTot;
	refine();
}

void HighsEquitable::refine(){
	if (refinements) findTarget();
	refinements++;
	double weight;
	colorsAdj.clear();
	colorsToSplit.clear();
	while (!S.empty()){
		int r = S.top();
		S.pop();
		SCheck[r] = false;
		bool var = r < nCols ? 1 : 0;
		for (int i = 0; i < Csize[r]; ++i){
			int v = C[r][i];
			if (var){
				for (int j = Astart[v]; j < Astart[v + 1]; ++j){
					int w = Aindex[j] + nCols;
					weight = AvaluePos[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]].push_back(w);
						Asize[color[w]]++;
					}
					if (!isAdj[color[w]]){
						isAdj[color[w]] = true;
						colorsAdj.push_back(color[w]);
					}
					if (cdeg[w] > maxcdeg[color[w]])
						maxcdeg[color[w]] = cdeg[w];			
				}
			}
			else{
				v -= nCols;
				for (int j = ARstart[v]; j < ARstart[v + 1]; ++j){
					int w = ARindex[j];
					weight = ARvaluePos[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]].push_back(w);
						Asize[color[w]]++;
					} 
					if (!isAdj[color[w]]){
						isAdj[color[w]] = true;
						colorsAdj.push_back(color[w]);
					}
					if (cdeg[w] > maxcdeg[color[w]])
						maxcdeg[color[w]] = cdeg[w];			
				}
			}
		}
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			if (Csize[c] != Asize[c])
				mincdeg[c] = 0;
			else{
				mincdeg[c] = maxcdeg[c];
				for (int j = 0; j < Asize[c]; ++j){
					int w = A[c][j];
					if (cdeg[w] < mincdeg[c]) mincdeg[c] = cdeg[w];
				}
			}
		}
		colorsToSplit.clear();
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			if (mincdeg[c] < maxcdeg[c])
				colorsToSplit.push_back(c);
		}
		sort(colorsToSplit.begin(), colorsToSplit.end());
		for (int i = 0; i < colorsToSplit.size(); ++i){
			int s = colorsToSplit[i];
			splitColor(s);
		}
		for (int i = 0; i < colorsAdj.size(); ++i){
			int c = colorsAdj[i];
			for (int j = 0; j < Asize[c]; ++j){
				int w = A[c][j];
				cdeg[w] = 0;
				numEdges[w] = 0;
			}
			maxcdeg[c] = 0;
			isAdj[c] = false;
			A[c].clear();
			Asize[c] = 0;
		}
		colorsAdj.clear();
	}			
	// cout << "\n\n Partition \n\n" << endl;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1) 
			isolates[C[i].front()] = true;
	}
	// for (int i = 0; i < coeff.size(); ++i){
	// 	cout << "node: " << i << ": ";
	// 	for (int j = 0; j < coeff[i].size(); ++j){
	// 		cout << coeff[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	packVectors();
}

void HighsEquitable::splitColor(int s){
	bool var = (s < nCols) ? true : false;
	set<double> cdegCopy;
	vector<int> colorFreq(nTot, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	degSumColor.insert(pair<double, int>(mincdeg[s], s));
	colorFreq[0] = Csize[s] - Asize[s];
	for (int i = 0; i < Asize[s]; ++i){
		int w = A[s][i];
		if (var){
			ret = degSumColor.insert(pair<double, int>(cdeg[w], vCol));
			if (ret.second){
				if (parentPartition[s] != -1) parentPartition[vCol] = parentPartition[s];
				else parentPartition[vCol] = s;
				vCol++;
			}
		}
		else{
			ret = degSumColor.insert(pair<double, int>(cdeg[w], cCol));
			cCol += ret.second == 1;
		}
		colorFreq[cdeg[w]]++;
	}
	int b = distance(colorFreq.begin(), max_element(colorFreq.begin(), colorFreq.end()));
	int instack = (SCheck[s]) ? 1 : 0;
	for(map<double, int>::iterator it = degSumColor.begin(); it != degSumColor.end(); ++it){
		//coeff[r].push_back(it->second);
		if (it->first == mincdeg[s]){
			if (!instack && it->first != b){
				S.push(it->second);
				SCheck[s] = true;
			}
		}
		else{
			if (it->first != b){
				S.push(it->second);
				SCheck[it->second] = true;
			}
		}
	}
	for (int i = 0; i < Asize[s]; ++i){
		int w = A[s][i];
		if (degSumColor[cdeg[w]] != s){
			C[s].erase(remove(C[s].begin(), C[s].end(), w), C[s].end());
			Csize[s]--;
			C[degSumColor[cdeg[w]]].push_back(w);
			Csize[degSumColor[cdeg[w]]]++;
			color[w] = degSumColor[cdeg[w]];
		}
	}
}

void HighsEquitable::handleNegatives(){
    vector<double>::iterator min = min_element(Avalue.begin(), Avalue.end());
    if (*min < 0){
        for (int i = 0; i < Avalue.size(); ++i){
            AvaluePos[i] = Avalue[i] + (-*min) + 1;
            ARvaluePos[i] = ARvalue[i] + (-*min) + 1;
        }
    }
    else{
        for (int i = 0; i < Avalue.size(); ++i){
            AvaluePos[i] = Avalue[i];
            ARvaluePos[i] = ARvalue[i];
        }
    }
}

void HighsEquitable::transpose(){
	int AcountX = Astart[nCols];
	ARindex.resize(AcountX);
	ARvalue.resize(AcountX);
	ARstart.resize(nRows + 1);
	vector<int> AR_Nend(nRows);
	for (int k = 0; k < AcountX; ++k) AR_Nend[Aindex[k]]++;
	for (int i = 1; i <= nRows; ++i) ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
	for (int i = 0; i < nRows; ++i) AR_Nend[i] = ARstart[i];
	for (int col = 0; col < nCols; ++col){
		for (int k = Astart[col]; k < Astart[col + 1]; ++k){
			int row = Aindex[k];
			int put = AR_Nend[row]++;
			ARindex[put] = col;
			ARvalue[put] = Avalue[k];
		}
	}
}

void HighsEquitable::findTarget(){
    for (int col = 0; col < nTot; ++col){
		if (Csize[col] > 1){
			// int isolated = i;
			isolate(col);
			return;
		}
	}
}

void HighsEquitable::isolate(int s){
	parentPartition.assign(nCols, -1);
	previousRowColoring.clear();
	previousColumnColoring.clear();
	prevC = C;
	for (int i = 0; i < nCols; ++i)
		previousColumnColoring.push_back(color[i]);
	for (int i = 0; i < nRows; ++i)
		previousRowColoring.push_back(color[i + nCols]);
	vector<int> temp1;
	vector<int> temp2;
	for (int i = 0; i < C[s].size(); ++i)
		if (!i) {temp1.push_back(C[s][i]); color[C[s][i]] = s;}
		else {temp2.push_back(C[s][i]); color[temp2[i - 1]] = vCol; Csize[vCol]++;}
	C[s] = temp1;
	C[vCol] = temp2;
	Csize[s] = 1;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1)
			isolates[C[i][0]] = true;
	}
	parentPartition[vCol] = s;
	SCheck[vCol] = true;
	S.push(vCol);
	vCol++;
}

bool HighsEquitable::isDiscrete(){
	for (int i = 0; i < nTot; ++i){
		if (Csize[i] > 1)
			return false;
	}
	return true;
}

void HighsEquitable::packVectors(){
	for (int i = 0; i < Aindex.size(); ++i)
		AindexP[i] = color[Aindex[i] + nCols];
}

// /* Sets labels for partition information (saucy style) */
// void HighsEquitable::setLabel(int index, int value){
// 	coloring.lab[index] = value;
// 	coloring.unlab[value] = index;
// }

// /* Add singletons and nonsingletons to induction list
// for refinement induction */
// void HighsEquitable::addInduce(int who){
// 	if (!coloring.clen[who]) {
// 		sinduce[nsinduce++] = who;
// 	}
// 	else {
// 		ninduce[nninduce++] = who;
// 	}
// 	indmark[who] = 1;
// }

// void HighsEquitable::fixFronts(int cf, int ff){
// 	int i, end = cf + coloring.clen[cf];
// 	for (i = ff; i <= end; ++i) {
// 		coloring.cfront[coloring.lab[i]] = cf;
// 	}
// }

// /* Allocate color storage */
// void HighsEquitable::colorAlloc(){
// 	coloring.lab = ints(nTot);
// 	coloring.unlab = ints(nTot);
// 	coloring.clen = ints(nTot);
// 	coloring.cfront = zeros(nTot);
// 	ccount = zeros(nTot);
// 	sinduce = ints(nTot);
// 	ninduce = ints(nTot);
// 	indmark = bits(nTot);
// 	count = ints(nTot);
// 	conncnts = ints(nTot);
// 	clist = ints(nTot);
// 	prevnon = ints(nTot + 1);
// 	nextnon = ints(nTot + 1) + 1;
// 	stuff = bits(nTot + 1);
// 	junk = ints(nTot);\
// 	bucket = ints(nTot + 2);
// }

// /* Going to try saucy refinement style */
// int HighsEquitable::refineSaucy(){
// 	int front;
// 	while(true){
// 		if (atTerminal()){
// 			clearRefine();
// 			return 1;
// 		}
// 		if (nsinduce){
// 			front = sinduce[--nsinduce];
// 			indmark[front] = 0;
// 			if (!refSingletonUndirected(front)) break;
// 		}
// 		else if (nninduce){
// 			front = ninduce[--nninduce];
// 			indmark[front] = 0;
// 			// if (!refNonsingleUndirected(front)) break;
// 		}
// 		else{
// 			return 1;
// 		}

// 	}
// }

// bool HighsEquitable::atTerminal(){
// 	return nsplits == nTot;
// }

// void HighsEquitable::clearRefine(){
// 	int i;
// 	for (i = 0; i < nninduce; ++i)
// 		indmark[ninduce[i]] = 0;
// 	for (i = 0; i < nsinduce; ++i)
// 		indmark[sinduce[i]] = 0;
// 	nninduce = nsinduce = 0;
// }

// void HighsEquitable::swapLabels(int a, int b){
// 	int tmp = coloring.lab[a];
// 	setLabel(a, coloring.lab[b]);
// 	setLabel(b, tmp);
// }

// void HighsEquitable::moveToBack(int k){
// 	int cf = coloring.cfront[k];
// 	int cb = cf + coloring.clen[cf];
// 	int offset = conncnts[cf]++;
// 	swapLabels(cb - offset, coloring.unlab[k]);
// 	if (!offset) clist[csize++] = cf;
// }

// void HighsEquitable::dataMark(int k){
// 	int cf = coloring.cfront[k];
// 	if (coloring.clen[cf]) moveToBack(k);
// }

// int HighsEquitable::refSingleton(int* adj, int* edg, int cf){
// 	int i, k = coloring.lab[cf];
// 	for (i = adj[k]; i < adj[k + 1]; ++i)
// 		dataMark(edg[i]);
// 		return refineCell();
// }

// int HighsEquitable::refNonsingle(int* adj, int*edg, int cf){
// 	int i, j, k, ret;
// 	int cb = cf + coloring.clen[cf];
// 	int size = cb - cf + 1;
// 	if (cf == cb)
// 		return refSingleton(adj, edg, cf);
// 	memcpy(junk, coloring.lab + cf, size * sizeof(int));
// 	for (i = 0; i < size; ++i){
// 		k = junk[i];
// 		for (j = adj[k]; j != adj[k + 1]; ++j)
// 			// dataCount(edg[j]); // TO DO: define dataCount function and edit it for lps
// 	}
// 	ret = refineCells();
// 	for (i = cf; i <= cb; ++i){
// 		k = coloring.lab[i];
// 		for (j = adj[k]; j != adj[k + 1]; ++j)
// 			ccount[edg[j]] = 0;
// 	}
// 	return ret;
// }

// void HighsEquitable::lp2Graph(){
// 	adj = ints(nTot + 1);
// 	edg = ints(2 * nnz);
// 	wt = doubles(2 * nnz);
// 	// Count adjacencies
// 	for (int i = 0; i < nCols; ++i){
// 		for (int j = Astart[i]; j < Astart[i + 1]; ++j){
// 			adj[i]++; adj[Aindex[j] + nCols]++; 
// 		}
// 	}
// 	// Insert adjacencies (sparse storage)
// 	int idx = 0;
// 	for (int i = 0; i < nCols; ++i){
// 		for (int j = Astart[i]; j < Astart[i + 1]; ++j){
// 			edg[j] = Aindex[j] + nCols;
// 			wt[j] = Avalue[j];
// 		}
// 	}
// 	for (int i = 0; i < nRows; ++i){
// 		for (int j = ARstart[i]; j < ARstart[i + 1]; ++j){
// 			edg[j + nnz] = ARindex[j];
// 			wt[j + nnz] = Avalue[j];
// 		}
// 	}
// 	fixAdj();
// }

// void HighsEquitable::fixAdj(){
// 	int i, val, sum;
// 	val = adj[0]; sum = 0; adj[0] = 0;
// 	for (i = 1; i <= nTot; ++i){
// 		sum += val;
// 		val = adj[i];
// 		adj[i] = sum;
// 	} 
// }

// int HighsEquitable::refineCell(){
// 	int i, cf, ret = 1;
// 	if (lev > 1) introsort(clist, csize); // need to define lev and introsort
// 	for (i = 0; ret && i < csize; ++i){
// 		cf = clist[i];
// 		ret = refSingleCell(cf);
// 	}
// 	for (i = 0; i < csize; ++i){
// 		cf = clist[i];
// 		conncnts[cf] = 0;
// 	}
// 	csize = 0;
// 	return ret;
// }

// int HighsEquitable::refineCells(){
// 	int i, cf, ret = 1;
// 	if (lev > 1) introsort(clist, csize); // need to define lev and introsort
// 	for (i = 0; ret && i < csize; ++i){
// 		cf = clist[i];
// 		ret = refNonsingleCell(cf);
// 	}
// 	for (i = 0; i < csize; ++i){
// 		cf = clist[i];
// 		conncnts[cf] = 0;
// 	}
// 	csize = 0;
// 	return ret;
// }

// int HighsEquitable::refSingleCell(int cf){
// 	int zcnt = coloring.clen[cf] + 1 - conncnts[cf];
// 	// return maybeSplit(cf, cf + zcnt); // TO DO define maybe split
// 	return 1;
// }

// int HighsEquitable::refNonsingleCell(int cf){
// 	int cnt, i, cb, nzf, ff, fb, bmin, bmax;
// 	cb = cf + coloring.clen[cf];
// 	nzf = cb - conncnts[cf] + 1;
// 	ff = nzf;
// 	cnt = ccount[coloring.lab[ff]];
// 	count[ff] = bmin = bmax = cnt;
// 	bucket[cnt] = 1;
// 	while (++ff <= cb){
// 		cnt = ccount[coloring.lab[ff]];
// 		while (bmin > cnt) bucket[--bmin] = 0;
// 		while (bmax < cnt) bucket[++bmax] = 0;
// 		++bucket[cnt];
// 		count[ff] = cnt;
// 	}
// 	if (bmin == bmax) return 1;
// 	ff = fb = nzf;
// 	for (i = bmin; i <= bmax; ++i, ff = fb){
// 		if (!bucket[i]) continue;
// 		fb = ff + bucket[i];
// 		bucket[i] = fb;
// 	}
// 	for (i = nzf; i <= cb; ++i)
// 		junk[--bucket[count[i]]] = coloring.lab[i];
// 	for (i = nzf; i <= cb; ++i)
// 		setLabel(i, junk[i]);
// 	for (i = bmax; i > bmin; --i){
// 		ff = bucket[i];
// 		if (ff && !split(cf, ff)) return 0; // TO DO: define split
// 	}
// 	return maybeSplit(cf, bucket[bmin]);
// }

// void HighsEquitable::introsort(int* a, int n){
// 	introsortLoop(a, n, 2 * logBase2(n));
// 	insertionSort(a, n);
// }

// void HighsEquitable::introsortLoop(int* a, int n, int lim){
// 	int p;
// 	while (n > 16){
// 		if (lim == 0){
// 			heapSort(a, n);
// 			return;
// 		}
// 		--lim;
// 		p = partition(a, n, median(a[0], a[n/2], a[n - 1]));
// 		introsortLoop(a + p, n - p, lim);
// 		n = p;
// 	}
// }

// int HighsEquitable::median(int a, int b, int c){
// 	if (a <= b) {
// 		if (b <= c) return b;
// 		if (a <= c) return c;
// 		return a;
// 	}
// 	else {
// 		if (a <= c) return a;
// 		if (b <= c) return c;
// 		return b;
// 	}
// }

// int HighsEquitable::partition(int* a, int n, int m){
// 	int f = 0, b = n;
// 	for (;;){
// 		while (a[f]) ++f;
// 		do --b; while (m <= a[b]);
// 		if (f < b){
// 			swap(a, f, b);
// 			++f;
// 		}
// 		else break;
// 	}
// 	return f;
// }

// int HighsEquitable::logBase2(int n){
// 	int k = 0;
// 	while (n > 1){
// 		++k;
// 		n >>= 1;
// 	}
// 	return k;
// }

// void HighsEquitable::swap(int* a, int x, int y){
// 	int tmp = a[x];
// 	a[x] = a[y];
// 	a[y] = tmp;
// }

// void HighsEquitable::heapSort(int* a, int n){
// 	int i;
// 	for (i = 1; i < n; ++i){
// 		siftUp(a - 1, i + 1);
// 	}
// 	--i;
// 	while(i > 0){
// 		swap(a, 0, i);
// 		siftDown(a - 1, i--);
// 	}
// }

// void HighsEquitable::siftUp(int* a, int k){
// 	int p;
// 	do{
// 		p = k/2;
// 		if(a[k] <= a[p]){
// 			return;
// 		}
// 		else{
// 			swap(a, k, p);
// 			k = p;
// 		}
// 	} while (k > 1);
// }

// void HighsEquitable::siftDown(int* a, int n){
// 	int p = 1, k = 2;
// 	while (k <= n){
// 		if (k < n && a[k] < a[k + 1]) ++k;
// 		if (a[p] < a[k]){
// 			swap(a, p, k);
// 			p = k;
// 			k = 2 * p;
// 		}
// 		else{
// 			return;
// 		}
// 	}
// }

// void HighsEquitable::insertionSort(int* a, int n){
// 	int i, j, k;
// 	for (i = 1; i < n; ++i){
// 		k = a[i];
// 		for (j = i; j > 0 && a[j - 1] > k; --j){
// 			a[j] = a[j - 1];
// 		}
// 		a[j] = k;
// 	}
// }

// int HighsEquitable::refSingletonUndirected(int cf){
// 	return refSingleton(adj, edg, cf);
// }

// int HighsEquitable::refNonsingleUndirected(int cf){
// 	return refNonsingle(adj, edg, cf);
// }

void HighsEquitable::amorph_print_automorphism(
    int n, const int *gamma, int nsupp, const int *support,
    struct amorph_graph *g, char *marks )
{
    int i, j, k;

    // We presume support is already sorted
    for (i = 0; i < nsupp; ++i) {
        k = support[i];

        // Skip elements already seen
        if (marks[k]) continue;

        // Start an orbit
        marks[k] = 1;
        printf( "(%s", g->var_names[k] );
        //printf("(%d", k);

        // Mark and notify elements in this orbit
        for (j = gamma[k]; j != k; j = gamma[j]) {
            marks[j] = 1;
            printf( " %s", g->var_names[j] );
            //printf(" %d", j);
        }

        // Finish off the orbit
        putchar(')');
    }
    putchar('\n');

    // Clean up after ourselves
    for (i = 0; i < nsupp; ++i) {
        marks[support[i]] = 0;
    }
}

int HighsEquitable::on_automorphism( int n, const int *gamma, int k, int *support, void *arg )
{
	// marks = (char*)calloc(nTotal, sizeof(char));
    // struct amorph_graph *g = (struct amorph_graph *)arg;
    // if( !quiet_mode )
    // {
    //     qsort_integers( support, k );
    //     //if( gap_mode )
    //     //{
    //     //    putchar( !first ? '[' : ',' );
    //     //    putchar( '\n' );
    //     //    first = 1;
    //     //}
    //     g->consumer( n, gamma, k, support, g, marks );
    // }
    return !timeout_flag;
}