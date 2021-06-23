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
	// Allocation for lpPartition class
	partition = (struct lpPartition *)calloc( (1), sizeof(struct lpPartition) );
	partition->cell.resize(nTot);
  	partition->cellFront.resize(nTot);
  	partition->cellSize.resize(nTot);
  	partition->labels.resize(nTot);
	partition->row.resize(nRows);
	partition->col.resize(nCols);
	partition->cellToCol.assign(nCols, -1);
	partition->cellToRow.assign(nRows, -1);
	partition->colsToReps.assign(nCols, -1);
	partition->repsToCols.assign(nCols, -1);
	partition->rowsToReps.assign(nRows, -1);
	partition->repsToRows.assign(nRows, -1);
	partition->parents.assign(nTot, -1);
}

int HighsEquitable::init_fixadj1( int n, int *adj ){
	int val, sum, i;
    /* Translate adj values to real locations */
    val = adj[0]; sum = 0; adj[0] = 0;
    for (i = 1; i < n; ++i) {
        sum += val;
        val = adj[i];
		adj[i] = sum;
    }
    return sum + val;
}

void HighsEquitable::init_fixadj2( int n, int e, int *adj ){
	int i;
	/* Translate again-broken sizes to adj values */
	for (i = n-1; i > 0; --i)
		adj[i] = adj[i-1];
	adj[0] = 0;
	adj[n] = e;
}

void HighsEquitable::amorph_print_automorphism(
    int n, const int *gamma, int nsupp, const int *support,
    struct amorph_graph *g, char *marks ){
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
    for (i = 0; i < nsupp; ++i)
        marks[support[i]] = 0;
}

int HighsEquitable::on_automorphism( int n, const int *gamma, int k, int *support, void *arg ){
    return !timeout_flag;
}

void HighsEquitable::lp2Graph(){
	int i, j, tempk, tempj, nColor = 0;
	int w, *aout, *ain, *eout, *ein, *colors, *wout, *win;
	double sum;
	// Use to color vertices by rhs, bounds, and obj coeff
	map<double, int> edgeColors;
	map<tuple<double, double>, int> rowColors;
	map<tuple<double, double, double>, int > colColors;
	aout = (int *)calloc( (nTot+1), sizeof(int) );
    eout = (int *)malloc( 2 * nnz * sizeof(int) );
    wout = (int *)malloc( 2 * nnz * sizeof(int) ); /* weight data */
    colors = (int *)malloc( nTot * sizeof(int) );
	char** var_names = (char**)malloc(nTot * sizeof(char*));
	//marks = (char*)calloc(nTot, sizeof(char));
	// Assign graph to temp pointers so it gets filled in
	g = (struct amorph_graph *)malloc(sizeof(struct amorph_graph));
	g->sg.nCols = nCols;
	g->sg.n = nTot;
	g->sg.e = nnz;
	g->sg.adj = aout;
	g->sg.edg = eout;
	g->sg.wght = wout;
	g->colors = colors;
	g->var_names = var_names;
	g->consumer = amorph_print_automorphism;
	// Connect in to out
	ain = aout;
	ein = eout;
	win = wout;
	// Fill initial colors
	pair<map<tuple<double, double, double>, int>::iterator, bool> ret;
	for (i = 0; i < nCols; ++i)
		if (colColors.insert(pair<tuple<double, double, double>, int>(
			make_tuple(colLower[i], colUpper[i], colCost[i]), nColor)).second)
			++nColor;
	for (i = 0; i < nRows; ++i)
		if (rowColors.insert(pair<tuple<double, double>, int>(
			make_tuple(rowLower[i], rowUpper[i]), nColor)).second)
			++nColor;
	nColor = 0;
	for (i = 0; i < nnz; ++i)
		if (edgeColors.insert(pair<double,int>(Avalue[i], nColor)).second)
			++nColor;
	g->sg.w = edgeColors.size();
	// Fill in colors
	for (i = nCols - 1; i >= 0; --i)
		colors[i] = colColors.find(make_tuple(colLower[i], colUpper[i], colCost[i]))->second;
	for (i = nRows - 1; i >= 0; --i)
		colors[i + nCols] = rowColors.find(make_tuple(rowLower[i], rowUpper[i]))->second;
	// Fill adj, should be the same as Astart and ARstart + nnz
	for(j = 0; j < nCols; ++j){
        for(i = Astart[j]; i < Astart[j+1]; ++i){
            ++ain[Aindex[i] + nCols]; ++aout[j];
        }
    }
    /* Fix that */
    init_fixadj1( nTot, aout );
	/* Insert adjacencies */
    for(j = 0; j < nCols; ++j){
        for(i = Astart[j]; i < Astart[j+1]; ++i){
            tempk = ain[Aindex[i] + nCols]++; tempj = aout[j]++;
            eout[tempj] = Aindex[i] + nCols; ein[tempk] = j;
            w = edgeColors.find(Avalue[i])->second;
            wout[tempj] = w; win[tempk] = w;
        }
    }
	// Fix edges and weights
	init_fixadj2(nTot, 2 * nnz, aout);	
	// Fill in var_names
	for (i = 0; i < nCols; ++i){
		tempk = colNames[i].length() + 1;
		var_names[i] = (char*)malloc(tempk*sizeof(char));
		strcpy(var_names[i], colNames[i].c_str());
	}
	for(i = 0; i < nRows; ++i){
        tempk = rowNames[i].length() + 1;
        var_names[i + nCols] = (char*)malloc(tempk * sizeof(char));
        strcpy(var_names[i + nCols], rowNames[i].c_str());
    }
}

void HighsEquitable::doSaucyEquitable(){
	int i, j, k = 0, maxSplit = 0;
	s = saucy_alloc(nTot, g->sg.w); // TO DO: add second argument to this function
	partitions = (struct eq_part *)calloc( (1), sizeof(struct eq_part) );
    // for( i = 0; i <+1; ++i)
    // {
    //     partitions[i].target = -2;
    //     partitions[i].labels = (int *)calloc(nTot, sizeof(int));
    //     partitions[i].fronts = (int *)calloc(nTot, sizeof(int));
	// 	partitions[i].parents = (int *)calloc(nTot, sizeof(int));
    // }
	partitions[0].target = -2;
	partitions[0].labels = (int *)calloc(nTot, sizeof(int));
	partitions[0].fronts = (int *)calloc(nTot, sizeof(int));
	partitions[0].parents = (int *)calloc(nTot, sizeof(int));
	saucy_search(s, &g->sg, 0, g->colors, on_automorphism, g, &stats, partitions); 
	int iter = stats.iter;
	// for (i = 0; i < iter; ++i)
	// 	if (partitions[i].nsplits > maxSplit) maxSplit = partitions[i].nsplits;
	// for (j = 0; j < i; ++j)
	// 	if (partitions[j].nsplits != maxSplit) partitions[j].nsplits = maxSplit;
	saucy_free(s);
}
// Refine function to call everything thing else
lpPartition* HighsEquitable::refine(){
	lp2Graph();
	doSaucyEquitable();
	convertToLpPartition();
	return partition;
	// return partitions;
} 

// Convert the saucy partition to lpPartition
void HighsEquitable::convertToLpPartition(){
	int i = 0, c = 0, r = 0, rep, cRep, pRep, v, pf, cf, colCounter = 0, rowCounter = 0;
	map<int, int> fCell;
	vector<bool> cellMarked(nTot, false);
	// map fronts to colors and count numCol_ and numRow_
	int numColCell = 0;
	for (i = 0; i < nTot; ++i){
		if (fCell.insert(pair<int, int>(partitions[0].fronts[i], c)).second){
		++c;
		++partition->numTot_;
		}
	}
	for (i = 0; i < nTot; ++i){
		partition->labels[i] = partitions[0].labels[i];
		if (i < nCols) partition->parents[i] = partitions[0].parents[i];
		if (partition->parents[i] > -1) partition->numResCol_++;
		partition->cell[i] = fCell.find(partitions[0].fronts[i])->second;
		partition->cellFront[fCell.find(partitions[0].fronts[i])->second] = partitions[0].fronts[i];
		++partition->cellSize[fCell.find(partitions[0].fronts[i])->second];
	}
	/* New mapping code: Go through col by col in ascending order,
	when a cell is counted it is marked as true, go to next col that
	that is in a new cell and mark that cell */
	for (i = 0; i < nTot; ++i){
		c = partition->cell[i];
		if (cellMarked[c]){ 
		i < nCols ? partition->col[i] = partition->cellToCol[c] : partition->row[i - nCols] = partition->cellToRow[c - partition->numCol_];
		continue;
		}
		if (i < nCols){
		partition->repsToCols[i] = partition->numCol_;
		partition->colsToReps[partition->numCol_] = i;
		partition->cellToCol[c] = partition->numCol_;
		partition->col[i] = partition->numCol_++;
		cellMarked[c] = true;
		}
		else{
		partition->repsToRows[i - nCols] = partition->numRow_;
		partition->rowsToReps[partition->numRow_] = i - nCols;
		partition->cellToRow[c - partition->numCol_] = partition->numRow_;
		partition->row[i - nCols] = partition->numRow_++;
		cellMarked[c] = true;
		}
	}
	partition->numTot_ = partition->numCol_ + partition->numRow_;
}

// Return number of refinements
int HighsEquitable::getNumRefinements(){
	int iter = stats.iter;
	return iter;
}
