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
	Xstart_.assign(numTot, 0);
	SCheck.assign(numTot, false);
	mincdeg.assign(numTot, 0);
	maxcdeg.assign(numTot, 0);
	cdeg.assign(numTot, 0);
	numEdges.assign(numTot, 0);
	isAdj.assign(numTot, false);
	isolates.assign(numTot, false);
	targets.assign(numTot, false);
	previousColumnColoring.assign(numCol, -1);
	previousRowColoring.assign(numRow, -1);
	partSize.assign(numCol, 0);
	previousPartSize.assign(numCol, 0);
	color.assign(numTot, 0);
	for (int i = 0; i < numTot; ++i){
		C.push_back(new list<int>);
		A.push_back(new forward_list<int>);
	}

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
		for (i = 0; i < Avalue.size(); ++i){
			AvalueCopy[i] = Avalue[i] + (-*min) + 1;
		}
	}
	else{
		for (i = 0; i < Avalue.size(); ++i)
			AvalueCopy[i] = Avalue[i];
	}
	// for (i = 0; i < numCol; ++i){
	// 	for (j = Astart[i]; j < Astart[i + 1]; ++j){
	// 		adjListLab[i].push_back(Aindex[j] + numCol);
	// 		adjListWeight[i].push_back(AvalueCopy[j]);
	// 		adjListWeightReal[i].push_back(Avalue[j]);
	// 		adjListLab[Aindex[j] + numCol].push_back(i);
	// 		adjListWeight[Aindex[j] + numCol].push_back(AvalueCopy[j]);
	// 		adjListWeightReal[Aindex[j] + numCol].push_back(Avalue[j]);
	// 	}
	// }
}

void HighsEquitable::initialRefinement(){
	int i;
	int numParts = 0;
	int varColor = 0;
	int conColor = numCol;
	set<double> Rhs;
	vector<double> Rhs_;
	set<pair<double, double> > rhs;
	set<tuple<double, double, double > > objBounds;
	// Prefill initial parts vector
	initialParts.assign(numTot, 0);

	// Prefill right hand side sets
	pair<set<pair<double, double> >::iterator, bool> retRow;
	for (i = 0; i < numRow; ++i){
		retRow = rhs.inset(pair<double, double>(rowLower[i], rowUpper[i]));
		if (retRow.second){
			initialParts[i + numCol] = conColor;
			conColor++;
		}
	}

	// Peace together var bounds and obj coefficients
	pair<set<tuple<double, double, double> >::iterator, bool> retCol;
	for (i = 0; i < numCol; ++i){
		retCol = objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i]));
		if (retCol.second){
			initialParts[i] = varColor;
			S.push(varColor);
			SCheck[varColor] = true;
			varColor++;
		}
	}

	// Define number of colors that have been used so far
	vCol = varColor;
	cCol = conColor;

	// Store initial refinement
	for (i = 0; i < numTot; ++i){
		C[initialParts[i]].push_back(i);
		color[i] = initialParts[i];
	}
}

void HighsEquitable::refine(){
	if (refinements) findTarget();
	refinements++;
	linkingPairs.clear();
	columnColorReps.clear();
	double weight;
	colorsAdj.clear();
	colorsToSplit.clear();
	numParts = 0;
	while (!S.empty()){
		r = S.top();
		S.pop();
		SCheck[r] = false;
		for (vPointer = C[r].begin(); vPointer != C[r].end(); ++vPointer){
			v = *vPointer;
			for (j = Astart[v]; j < Astart[v + 1]; ++j){
				w = Aindex[j];
				weight = Avalue[j];
				cdeg[w] += weight;
				numEdges[w]++;
				if (numEdges[w] == 1)
					A[color[w]].push_back(w);
				if (!isAdj[colow[w]]){
					isAdj[color[w]] = true;
					colorsAdj.push_back(color[w]);
				}
				if (cdeg[w] > maxcdeg[color[w]])
					maxcdeg[color[w]] = cdeg[w];			
			}
		}
		for (cPointer = colorsAdj.begin(); cPointer!= colorsAdj.end(); ++cPointer){
			c = *cPointer;
			if (C[c].size() != A[c].size())
				mincdeg[c] = 0;
			else{
				mincdeg[c] = maxcdeg[c];
				for (wPointer = A[c].begin(); wPointer != A[c].end(); ++wPoiner){
					w = *wPointer;
					if (cdeg[w] < mincdeg[c]) mincdeg[c] = cdeg[w];
				}
			}
		}
		colorsToSplit.clear();
		for (cPointer = colorsAdj.begin(); cPointer != colorsAdj.end; ++cPointer){
			c = *cPointer;
			if (mincdeg[c] < maxcdeg[c])
				colorsToSplit.push_back(c);
		}
		colorsToSplit.sort();
		for (sPointer = colorsToSplit.begin(); sPointer != colorsToSplit.end(); ++sPointer){
			s = *sPointer;
			splitColor(s);
		}
		for (cPointer = colorsAdj.begin(); cPointer != colorsAdj.end(); ++cPointer){
			c = *cPointer;
			for (wPointer = A[c].begin(); wPointer != A[c].end(); ++wPointer){
				w = *wPointer;
				cdeg[w] = 0;
				numEdges[w] = 0;
			}
			maxcdeg[c] = 0;
			A[c].clear();
		}
		colorsAdj.clear();
	}			
	//cout << "\n\n Partition \n\n" << endl;
	for (int i = 0; i < C.size(); ++i){
		if (C[i].size() == 1) 
			isolates[C[i].front()] = true;
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
	bool varOrCon = (s < numCol) ? true : false;
	set<double> cdegCopy;
	vector<int> colorFreq(numTot, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	degSumColor.insert(pair<double, int>(mincdeg[s], s));
	colorFreq[0] = C[s].size() - A[s].size();
	for (wPointer = A[s].begin(); wPointer != A[s].end(); ++wPointer){
		w = *wPointer;
		if (varOrCon){
			ret = degSumColor.insert(pair<double, int>(cdeg[w], vCol));
			vCol += (ret.second == 1);
		}
		else{
			ret = degSumColor.insert(pair<double, int>(cdeg[w], cCol));
			cCol += (ret.second == 1);
		}
		colorFreq[degSumColor[cdeg[w]]]++;
	}
	if (varOrCon){
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
	}
	for (wPointer = A[s].begin(); wPointer != A[s].end(); ++wPointer){
		w = *wPointer;
		if (degSumColor[cdeg[w]] != s){
			C[s].remove(w);
			C[degSumColor[cdeg[w]]].push_back(w);
			color[w] = degSumColor[cdeg[w]];
		}
	}
}

void HighsEquitable::findTarget(){
	for (int i = 0; i < numTot; ++i){
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
	C[color[i]].erase(remove(C[oldCol].begin(), C[oldCol].end(), i), C[oldCol].end());
	C[newCol].push_back(i);
	color[i] = newCol;
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
	for (int i = 0; i < numCol; ++i){
		if (!isolates[i]) return false;
	}
	return true;
}