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
	coeff.resize(numTot);
	Xstart_.assign(Astart.begin(), Astart.end());
	Xvalue_.assign(Avalue.size(), 0);
	Xindex_.assign(Aindex.size(), 0);
	AvalueCopy.assign(Avalue.size(), 0);
	ARvalueCopy.assign(Avalue.size(), 0);
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
	Csize.assign(numTot, 0);
	Asize.assign(numTot, 0);
	for (int i = 0; i < numTot; ++i){
		C.push_back(new list<int>());
		A.push_back(new forward_list<int>());
	}

	// Initial refinement
	handleNegatives();
	createRowCopy();
	initialRefinement();
	refine();
}

void HighsEquitable::handleNegatives(){
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
}

void HighsEquitable::createRowCopy(){
    int AcountXSub = Astart[numCol];
    ARindex.resize(AcountXSub);
    ARvalueCopy.resize(AcountXSub);
    // Build row copy - pointers
    ARstart.assign(numRow + 1, 0);
    AR_Nend.assign(numRow, 0);
    for (int k = 0; k < AcountXSub; ++k)
        AR_Nend[Aindex[k]]++;
    for (int i = 1; i <= numRow; ++i)
        ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
    for (int i = 0; i < numRow; ++i)
        AR_Nend[i] = ARstart[i];
    // Build row copy - elements
    for (int iCol = 0; iCol < numCol; ++iCol) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; ++k) {
            int iRow = Aindex[k];
            int iPut = AR_Nend[iRow]++;
            ARindex[iPut] = iCol;
            ARvalueCopy[iPut] = AvalueCopy[k];
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
	set<pair<double, double> > rhs;
	set<tuple<double, double, double > > objBounds;
	// Prefill initial parts vector
	initialParts.assign(numTot, 0);

	// Prefill right hand side sets
	pair<set<pair<double, double> >::iterator, bool> retRow;
	rhs.insert(pair<double, double>(rowLower[0], rowUpper[0]));
	initialParts[numCol] = conColor;
	S.push(conColor);
	SCheck[conColor] = true;
	for (i = 1; i < numRow; ++i){
		retRow = rhs.insert(pair<double, double>(rowLower[i], rowUpper[i]));
		if (retRow.second){
			initialParts[i + numCol] = ++conColor;
			S.push(conColor);
			SCheck[conColor] = true;
		}
		else
			initialParts[i + numCol] = conColor;
	}
	conColor++;

	// Peace together var bounds and obj coefficients
	pair<set<tuple<double, double, double> >::iterator, bool> retCol;
	objBounds.insert(make_tuple(colCost[0], colLower[0], colUpper[0]));
	initialParts[0] = varColor;
	S.push(varColor);
	SCheck[varColor] = true;
	for (i = 1; i < numCol; ++i){
		retCol = objBounds.insert(make_tuple(colCost[i], colLower[i], colUpper[i]));
		if (retCol.second){
			initialParts[i] = ++varColor;
			S.push(varColor);
			SCheck[varColor] = true;
		}
		else
			initialParts[i] = varColor;
	}
	varColor++;

	// Define number of colors that have been used so far
	vCol = varColor;
	cCol = conColor;

	// Store initial refinement
	for (i = 0; i < numTot; ++i){
		C[initialParts[i]]->push_back(i);
		Csize[initialParts[i]]++;
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
		bool varOrCon = r < numCol ? 1 : 0;
		for (vPointer = C[r]->begin(); vPointer != C[r]->end(); ++vPointer){
			v = *vPointer;
			if (varOrCon){
				for (int j = Astart[v]; j < Astart[v + 1]; ++j){
					w = Aindex[j] + numCol;
					weight = AvalueCopy[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]]->push_front(w);
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
				v -= numCol;
				for (int j = ARstart[v]; j < ARstart[v + 1]; ++j){
					w = ARindex[j];
					weight = ARvalueCopy[j];
					cdeg[w] += weight;
					numEdges[w]++;
					if (numEdges[w] == 1){
						A[color[w]]->push_front(w);
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
		for (cPointer = colorsAdj.begin(); cPointer!= colorsAdj.end(); ++cPointer){
			c = *cPointer;
			if (Csize[c] != Asize[c])
				mincdeg[c] = 0;
			else{
				mincdeg[c] = maxcdeg[c];
				for (wPointer = A[c]->begin(); wPointer != A[c]->end(); ++wPointer){
					w = *wPointer;
					if (cdeg[w] < mincdeg[c]) mincdeg[c] = cdeg[w];
				}
			}
		}
		colorsToSplit.clear();
		for (cPointer = colorsAdj.begin(); cPointer != colorsAdj.end(); ++cPointer){
			c = *cPointer;
			if (mincdeg[c] < maxcdeg[c])
				colorsToSplit.push_front(c);
		}
		colorsToSplit.sort();
		for (sPointer = colorsToSplit.begin(); sPointer != colorsToSplit.end(); ++sPointer){
			s = *sPointer;
			splitColor(s);
		}
		for (cPointer = colorsAdj.begin(); cPointer != colorsAdj.end(); ++cPointer){
			c = *cPointer;
			for (wPointer = A[c]->begin(); wPointer != A[c]->end(); ++wPointer){
				w = *wPointer;
				cdeg[w] = 0;
				numEdges[w] = 0;
			}
			maxcdeg[c] = 0;
			isAdj[c] = false;
			A[c]->clear();
			Asize[c] = 0;
		}
		colorsAdj.clear();
	}			
	// cout << "\n\n Partition \n\n" << endl;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1) 
			isolates[C[i]->front()] = true;
		if (Csize[i]){
			// if (i < numCol) partSize[i] = C[i].size();
			// //cout << "color: " << i << endl;
			rep = C[i]->front();
			if (i < numCol && !targets[rep])
				columnColorReps.push_back(rep);
		}
	}
	// for (int i = 0; i < coeff.size(); ++i){
	// 	cout << "node: " << i << ": ";
	// 	for (int j = 0; j < coeff[i].size(); ++j){
	// 		cout << coeff[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	packVectors();
	collectLinkingPairs();
}

void HighsEquitable::splitColor(int s){
	bool varOrCon = (s < numCol) ? true : false;
	set<double> cdegCopy;
	vector<int> colorFreq(numTot, 0);
	map<double, int> degSumColor;
	pair<map<double, int>::iterator, bool> ret;
	degSumColor.insert(pair<double, int>(mincdeg[s], s));
	colorFreq[0] = Csize[s] - Asize[s];
	for (wPointer = A[s]->begin(); wPointer != A[s]->end(); ++wPointer){
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
	int b = distance(colorFreq.begin(), max_element(colorFreq.begin(), colorFreq.end()));
	int instack = (SCheck[s]) ? 1 : 0;
	for(map<double, int>::iterator it = degSumColor.begin(); it != degSumColor.end(); ++it){
		coeff[r].push_back(it->second);
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
	for (wPointer = A[s]->begin(); wPointer != A[s]->end(); ++wPointer){
		w = *wPointer;
		if (degSumColor[cdeg[w]] != s){
			C[s]->remove(w);
			Csize[s]--;
			C[degSumColor[cdeg[w]]]->push_back(w);
			Csize[degSumColor[cdeg[w]]]++;
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
	C[oldCol]->remove(i);
	Csize[oldCol]--;
	C[newCol]->push_back(i);
	Csize[newCol]++;
	color[i] = newCol;
	vCol++;
	for (int i = 0; i < C.size(); ++i){
		if (Csize[i] == 1)
			isolates[C[i]->front()] = true;
	}
	SCheck[oldCol] = true;
	S.push(oldCol);
}

void HighsEquitable::packVectors(){
	for (int i = 0; i < Aindex.size(); ++i){
		Xindex_[i] = color[Aindex[i] + numCol];
	}
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