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

void HighsEquitable::setup(const HighsLp& lp){
	// From the lp 
	numRow = (lp.numRow_);
	numCol = (lp.numCol_);
	nnz = (lp.nnz_);
	rowLower = (lp.rowLower_);
	rowUpper = (lp.rowUpper_);
	colUpper = (lp.colUpper_);
	colLower = (lp.colLower_);
	colCost = (lp.colCost_);
	Avalue = (lp.Avalue_);
	Aindex = (lp.Aindex_);
	Astart = (lp.Astart_);
	model_name = lp.model_name_;
	lp_name = lp.lp_name_;

	// Associated with equitable partition
	numTot = numCol + numRow;
	SCheck.assign(numTot, false);
	mincdeg.assign(numTot, 0);
	maxcdeg.assign(numTot, 0);
	cdeg.assign(numTot, 0);
	isAdj.assign(numTot, 0);
	isolates.assign(numTot, false);
	targets.assign(numTot, false);
	previousColumnColoring.assign(numCol, -1);
	previousRowColoring.assign(numRow, -1);
	partSize.assign(numCol, 0);
	previousPartSize.assign(numCol, 0);
	adjListLab.assign(numTot, vector<int>(0));
	adjListWeight.assign(numTot, vector<double>(0));
	C.assign(numTot, vector<int>(0));
	A.assign(numTot, vector<int>(0));

	// Initial refinement
	lp2Graph();
	initialRefinement();
	refine();
}

void HighsEquitable::lp2Graph(){
	vector<double> AvalueCopy(Avalue.size(), 0);
	vector<double>::iterator min = min_element(Avalue.begin(), Avalue.end());
	int i, j, k;
	if (*min < 0){
		for (i = 0; i < Avalue.size(); ++i)
			AvalueCopy[i] = Avalue[i] + -2*(*min);
	}
	else{
		for (i = 0; i < Avalue.size(); ++i)
			AvalueCopy[i] = Avalue[i];
	}
	for (i = 0; i < numCol; ++i){
		for (j = Astart[i]; j < Astart[i + 1]; ++j){
			adjListLab[i].push_back(Aindex[j] + numCol);
			adjListWeight[i].push_back(AvalueCopy[j]);
			adjListLab[Aindex[j] + numCol].push_back(i);
			adjListWeight[Aindex[j] + numCol].push_back(AvalueCopy[j]);
		}
	}
}

void HighsEquitable::initialRefinement(){
	int i;
	int numParts = 0;
	int varColor = 0;
	int conColor = numCol;
	set<double> Rhs;
	vector<double> Rhs_;
	set<double>::iterator rhsIdx;
	set<tuple<double, double, double > > objBounds;
	set<tuple<double, double, double > >::iterator objBoundsIdx;
	// Prefill initial parts vector
	for (i = 0; i < numTot; ++i)
		initialParts.push_back(0);

	// Prefill right hand side sets
	double rhs = 0;
	double absUpper = 0;
	double absLower = 0;
	for (i = 0; i < numRow; ++i){
		absUpper = fabs(rowUpper[i]);
		absLower = fabs(rowLower[i]);
		rhs = min(absUpper, absLower);
		Rhs_.push_back(rhs);
		Rhs.insert(rhs);
	}

	// Peace together var bounds and obj coefficients
	for (i = 0; i < numCol; ++i)
		objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i]));
	
	// Intial variable refinement
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

	// Initial constraint refinement
	for (rhsIdx = Rhs.begin(); rhsIdx != Rhs.end(); ++rhsIdx){
		for (int i = 0; i < numRow; ++i){
			if (Rhs_[i] == *rhsIdx){
				initialParts[i + numCol] = conColor;
			}
		}
		S.push(conColor);
		SCheck[conColor] = true;
		conColor++;
	}

	// Define number of colors that have been used so far
	vCol = varColor;
	cCol = conColor;

	// Store initial refinement
	for (i = 0; i < numTot; ++i){
		C[initialParts[i]].push_back(i);
		color.push_back(initialParts[i]);
	}
}

void HighsEquitable::refine(){
	if (refinements) findTarget();
	refinements++;
	linkingPairs.clear();
	columnColorReps.clear();
	int i, j, k, u, v, w, adjColor, rep;
	double weight;
	colorsAdj.clear();
	colorsToSplit.clear();
	numParts = 0;
	while (!S.empty()){
		r = S.top();
		S.pop();
		SCheck[r] = false;
		for (i = 0; i < C[r].size(); ++i){
			v = C[r][i];
			for (j = 0; j < adjListLab[v].size(); ++j){
				w = adjListLab[v][j];
				weight = adjListWeight[v][j];
				cdeg[w] += weight;
				isAdj[w]++;
				if (isAdj[w] == 1)
					A[color[w]].push_back(w);
				if (!(find(colorsAdj.begin(), colorsAdj.end(), color[w]) != colorsAdj.end()))
					colorsAdj.push_back(color[w]);
				if (cdeg[w] > maxcdeg[color[w]])
					maxcdeg[color[w]] = cdeg[w];			
			}
		}
		for (i = 0; i < colorsAdj.size(); ++i){
			adjColor = colorsAdj[i];
			if (C[adjColor].size() != A[adjColor].size())
				mincdeg[adjColor] = 0;
			else{
				mincdeg[adjColor] = maxcdeg[adjColor];
				for (j = 0; j < A[adjColor].size(); ++j){
					u = A[adjColor][j];
					if (cdeg[u] < mincdeg[adjColor]) mincdeg[adjColor] = cdeg[u];
				}
			}
		}
		colorsToSplit.clear();
		for (i = 0; i < colorsAdj.size(); ++i){
			adjColor = colorsAdj[i];
			if (mincdeg[adjColor] < maxcdeg[adjColor])
				colorsToSplit.push_back(adjColor);
		}
		sort(colorsToSplit.begin(), colorsToSplit.end());
		for (i = 0; i < colorsToSplit.size(); ++i){
			adjColor = colorsToSplit[i];
			splitColor(adjColor);
		}
		for (i = 0; i < colorsAdj.size(); ++i){
			u = colorsAdj[i];
			for (j = 0; j < A[u].size(); ++j){
				v = A[u][j];
				cdeg[v] = 0;
				isAdj[v] = 0;
			}
			maxcdeg[u] = 0;
			A[u].clear();
		}
		colorsAdj.clear();
		for (i = 0; i < C.size(); ++i)
			if (C[i].size() == 1) isolates[C[i].front()] = true;
	}
	//cout << "\n\n Partition \n\n" << endl;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size()){
			if (i < numCol) partSize[i] = C[i].size();
			//cout << "color: " << i << endl;
			rep = C[i].front();
			if (i < numCol && !targets[rep])
				columnColorReps.push_back(rep);
		}
	}
	collectLinkingPairs();
}

void HighsEquitable::splitColor(int s){
	int i,j,k,u,w,v;
	bool varOrCon = (s < numCol) ? true : false;
	set<double> cdegCopy;
	vector<int> colorFreq(numTot, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	degSumColor.insert(pair<double, int>(mincdeg[s], s));
	colorFreq[0] = C[s].size() - A[s].size();
	for (i = 0; i < A[s].size(); ++i){
		u = A[s][i];
		if (varOrCon){
			ret = degSumColor.insert(pair<double, int>(cdeg[u], vCol));
			vCol += (ret.second == 1);
		}
		else{
			ret = degSumColor.insert(pair<double, int>(cdeg[u], cCol));
			cCol += (ret.second == 1);
		}
	}
	for (i = 0; i < A[s].size(); ++i){
		u = A[s][i];
		colorFreq[degSumColor[cdeg[u]]]++;
	}
	int b = distance(colorFreq.begin(), max_element(colorFreq.begin(), colorFreq.end()));
	int instack = (SCheck[s]) ? 1 : 0;
	for(map<double, int>::iterator it = degSumColor.begin(); it != degSumColor.end(); ++it){
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
	for (i = 0; i < A[s].size(); ++i){
		u = A[s][i];
		if (degSumColor[cdeg[u]] != s){
			C[s].erase(remove(C[s].begin(), C[s].end(), u), C[s].end());
			C[degSumColor[cdeg[u]]].push_back(u);
			color[u] = degSumColor[cdeg[u]];
		}
	}
}

void HighsEquitable::findTarget(){
	for (int i = 0; i < numCol; ++i){
		if (!isolates[i]){
			isolated = i;
			targets[i] = true;
			isolate(i);
			return;
		}
	}
}

void HighsEquitable::isolate(int i){
	for (int j = 0; j < numCol; ++j)
		previousColumnColoring[j] = color[j];
	for (int j = numCol; j < numTot; ++j)
		previousRowColoring[j - numCol] = color[j];
	previousPartSize = partSize;
	prevC = C;
	// C[color[i]].erase(remove(C[color[i]].begin(), C[color[i]].end(), i), C[color[i]].end());
	int newCol = vCol;
	int oldCol = color[i];
	for (int j = 0; j < C[oldCol].size();){
		int var = C[oldCol][j];
		if (var != i){
			C[oldCol].erase(C[oldCol].begin() + j);
			C[newCol].push_back(var);
			color[var] = newCol;
		}
		else
			++j;
	}
	vCol++;
	for (int i = 0; i < C.size(); ++i){
		if (C[i].size() == 1)
			isolates[C[i].front()] = true;
	}
	SCheck[oldCol] = true;
	S.push(oldCol);
}

void HighsEquitable::collectLinkingPairs(){
	if (isolated == -1) return;
	int i, linkCnt = 0;
	int rep;
	int previousColorRep;
	int targ = isolated;
	int previousColorTarg = previousColumnColoring[targ];
	vector<int> notLinkedToIsolated;
	vector<int> temp;
	commonLinkers.clear();
	for (i = 0; i < columnColorReps.size(); ++i){
		rep = columnColorReps[i];
		previousColorRep = previousColumnColoring[rep];
		if (previousColorRep == previousColorTarg){
			linkingPairs.push_back(pair<int, int>(color[targ], color[rep]));
			temp.push_back(linkCnt);
			linkCnt++;
		}
		else{
			notLinkedToIsolated.push_back(rep);
		}
	}
	commonLinkers.insert(pair<int, vector<int> >(color[targ], temp));
	temp.clear();
	reverse(notLinkedToIsolated.begin(), notLinkedToIsolated.end());
	while (!notLinkedToIsolated.empty()){
		targ = notLinkedToIsolated.back();
		notLinkedToIsolated.pop_back();
		previousColorTarg = previousColumnColoring[targ];
		for (i = 0; i < notLinkedToIsolated.size();){
			rep = notLinkedToIsolated[i];
			previousColorRep = previousColumnColoring[rep];
			if (previousColorRep == previousColorTarg){
				linkingPairs.push_back(pair<int, int>(color[targ], color[rep]));
				temp.push_back(linkCnt);
				linkCnt++;
				notLinkedToIsolated.erase(notLinkedToIsolated.begin() + i);
			}
			else
				++i;
		}
		if (temp.size()){
			commonLinkers.insert(pair<int, vector<int> >(color[targ], temp));
			temp.clear();
		}
	}
} 

bool HighsEquitable::isPartitionDiscrete(){
	for (int i = 0; i < numTot; ++i){
		if (!isolates[i]) return false;
	}
	return true;
}