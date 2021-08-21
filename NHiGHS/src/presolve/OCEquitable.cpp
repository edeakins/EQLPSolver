#include "OCEquitable.h"
using namespace std;

void HighsOCEquitablePartition::runToDiscrete(){
    while (!discrete()) isolate();
}

void HighsOCEquitablePartition::isolate(){
    int targ = nextNon[0];
    int min = targ;
    int back = targ + partition->setLen[targ];
    // Additions, hopefully better
    front[targ + partition->setLen[targ]] = partition->setLen[targ];
    swapLabels(min, back);
    split(targ, back);
    refine();
}

void HighsOCEquitablePartition::refine(){
    int front;
    while (!S.empty()){
        if (discrete()){ 
            clear();
            break;
        }
        front = S.top(), S.pop();
        if (partition->setLen[front]){
            nIndSize--;
            indMark[front] = false;
            if (!refineNonsingleCell(front)) break;
        }
        else{ 
            sIndSize--;
            indMark[front] = false;
            if (!refineSinlgeCell(front)) break;
        }
    }
    clear();
}

// Refining cell was nonsingleton
bool HighsOCEquitablePartition::refineNonsingleCell(int sf){
    int i, j, k, edg, wght, sb, cf, ssize, e, w, 
    wcount = 0, ecount = 0, numAdj = 0, idx;
    bool ret = true;
    // sf = partition->front[node];
    sb = sf + partition->setLen[sf];
    ssize = sb - sf + 1;
    // Compute degree sums for all vertices adjacent to refining cell
    for (i = sf; i < sb + 1; ++i){
        j = partition->label[i];
        for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
            edg = g->edg[k]; wght = g->wght[k]; cf = partition->front[edg];
            if (!wghtFreq[wght]++) wcount++;
            vEdgColor[ecount] = wght;
            vEdg[ecount++] = edg;
        }
    }
    // Populate sparse start array for edge weights
    for (i = 0; i < wcount; ++i){
        wghtStart[i + 1] += wghtFreq[i];
        wghtOffset[i] = wghtStart[i];
    }
    // Populate sparse weight index array.  
    // Don't need sparse val array as the indices are stored
    // according to their weight.
    for (i = 0; i < ecount; ++i){
        w = vEdgColor[i];
        e = vEdg[i];
        idx = wghtOffset[w]++;
        wghtEdg[idx] = e;
    }

    for (i = 0; i < wcount; ++i){
        for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
            dataCount(wghtEdg[j]);
        }
        ret = ret && refineNonSingles();
        for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
            scount[wghtEdg[j]] = 0;
            wghtEdg[j] = 0;
            vEdgColor[j] = 0;
            vEdg[j] = 0;
        }
        wghtStart[i] = 0;
        wghtOffset[i] = 0;
        wghtFreq[i] = 0;
    }
    wghtStart[wcount] = 0;
    return ret;
}

// Refining cell is a singleton
bool HighsOCEquitablePartition::refineSinlgeCell(int sf){
    int i, j = partition->label[sf], k, edg, wght, sb, cf, 
    ssize, e, w, wcount = 0, ecount = 0, numAdj = 0, idx, cdeg;
    bool ret = true;
    maxDeg = 0; adjCellCnt = 0;
    for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
        edg = g->edg[k]; wght = g->wght[k]; cf = partition->front[edg];
        if (!wghtFreq[wght]++) wcount++;
        vEdgColor[ecount] = wght;
        vEdg[ecount++] = edg;
        // cDeg[edg] += wght + 1;
        if ((cDeg[edg] += wght + 1) > maxDeg) maxDeg = cDeg[edg];
        if (!nodeAdj[edg]++) refSize[cf]++;
        if (!cellAdj[cf]++) markCell(cf);
    }

    for (k = 0; k < adjCellCnt; ++k){
        splitCell(adjCell[k]);
        adjCell[k] = 0;
    }
    for (k = g->adj[j]; k < g->adj[j + 1]; ++k)
        if (cDeg[g->edg[k]]) cDeg[g->edg[k]] = 0;
    adjCellCnt = 0;
    // Populate sparse start array for edge weights
    for (i = 0; i < wcount; ++i){
        wghtStart[i + 1] += wghtFreq[i];
        wghtOffset[i] = wghtStart[i];
    }

    // Populate sparse weight index array.  
    // Don't need sparse val array as the indices are stored
    // according to their weight.
    for (i = 0; i < ecount; ++i){
        w = vEdgColor[i];
        e = vEdg[i];
        idx = wghtOffset[w]++;
        wghtEdg[idx] = e;
    }

    for (i = 0; i < wcount; ++i){
        for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
            dataMark(wghtEdg[j]);
        }
        ret = ret && refineSingles();
        for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
            scount[wghtEdg[j]] = 0;
            wghtEdg[j] = 0;
            vEdgColor[j] = 0;
            vEdg[j] = 0;
        }
        wghtStart[i] = 0;
        wghtOffset[i] = 0;
        wghtFreq[i] = 0;
    }
    wghtStart[wcount] = 0;
    return ret;
}

// Calls splitnonsingle
bool HighsOCEquitablePartition::refineNonSingles(){
    int i, cf, ret = true;
    // iterate over cells marked for refinement
    for (i = 0; ret && i < sListSize; ++i){
        cf = sList[i];
        ret = splitNonsingleCell(cf);
    }
    // Clear out connections
    for (i = 0; i < sListSize; ++i){
        cf = sList[i];
        conncnts[cf] = 0;
    }
    sListSize = 0;
    return ret;
}

// Calls split single
bool HighsOCEquitablePartition::refineSingles(){
    int i, cf, ret = true;
    // iterate over cells marked for refinement
    for (i = 0; ret && i < sListSize; ++i){
        cf = sList[i];
        ret = splitSinlgeCell(cf);
    }
    // Clear out connections
    for (i = 0; i < sListSize; ++i){
        cf = sList[i];
        conncnts[cf] = 0;
    }
    sListSize = 0;
    return ret;
}

// Split cell when refining cell was nonsingleton
bool HighsOCEquitablePartition::splitNonsingleCell(int cf){
    int cnt, i, cb, nzf, ff, fb, sMin, sMax;
    cb = cf + partition->setLen[cf];
    nzf = cb - conncnts[cf] + 1;
    ff = nzf;
    cnt = scount[partition->label[ff]];
    count[ff] = sMin = sMax = cnt;
    newSetSize[cnt] = 1;
    while (++ff <= cb){
        cnt = scount[partition->label[ff]];
        while (sMin > cnt) newSetSize[--sMin] = 0;
        while (sMax < cnt) newSetSize[++sMax] = 0;
        ++newSetSize[cnt];
        count[ff] = cnt;
    }
    if (sMin == sMax && cf == nzf) return true;
    ff = fb = nzf;
    for (i = sMin; i <= sMax; ++i, ff = fb){
        if (!newSetSize[i]) continue;
        fb = ff + newSetSize[i];
        newSetSize[i] = fb;
    }
    for (i = nzf; i <= cb; ++i)
        splitSet[--newSetSize[count[i]]] = partition->label[i];
    for (i = nzf; i <= cb; ++i)
        setLabel(i, splitSet[i]);
    for (i = sMax; i > sMin; --i){
        ff = newSetSize[i];
        if (ff && !split(cf, ff)) return false;
    }
    return possiblySplit(cf, newSetSize[sMin]);
}

// Split cell when refining cell was singleton
bool HighsOCEquitablePartition::splitSinlgeCell(int cf){
    int zcnt = partition->setLen[cf] + 1 - conncnts[cf];
    return possiblySplit(cf, cf + zcnt);
}

// Does the acutal splitting of a cell once new cell sizes are known
bool HighsOCEquitablePartition::split(int cf, int ff){
    // Track number of splits
    ++nSplits;
    // Do splitting
    int cb, fb;
    fb = ff - 1;
    cb = cf + partition->setLen[cf];
    partition->setLen[cf] = fb - cf;
    partition->setLen[ff] = cb - ff;
    // Fix the fronts
    fixFronts(ff, ff);
    // Add inducers
    if (indMark[cf] || partition->setLen[ff] < partition->setLen[cf])
        addInduce(ff);
    else   
        addInduce(cf);
    // Keep track of non singletons for new targets for inducing another level of refinement
    if (partition->setLen[ff]){
        prevNon[nextNon[cf + 1]] = ff;
        nextNon[ff + 1] = nextNon[cf + 1];
        prevNon[ff] = cf;
        prevNon[cf] = ff;
    }
    if (!partition->setLen[cf]){
        nextNon[prevNon[cf] + 1] = nextNon[cf + 1];
        prevNon[nextNon[cf + 1]] = prevNon[cf];
    }
    return true;
}

// Possibly does actual splitting when there is a 0 wght frequency
bool HighsOCEquitablePartition::possiblySplit(int cf, int ff){
    return cf == ff ? true : split(cf, ff);
}

void HighsOCEquitablePartition::allocatePartition(){
    // Convert LP to sparse graph rep
    lp2Graph();
    int i, j, max = 0, maxSetSize = 0, idx = 0;
    // maxSetSize = g->numCol_ > g->numRow_ ? g->numCol_ : g->numRow_;
    // Allocate partition storage
    partition = (struct OCpartition*)calloc(1, sizeof(struct OCpartition));
    partition->target = -1;
    partition->nplits = 0;
    partition->front.resize(g->numTot_);
    partition->label.resize(g->numTot_);
    partition->unlabel.resize(g->numTot_);
    partition->parent.resize(g->numTot_);
    partition->set.resize(g->numTot_);
    partition->setSize.resize(g->numTot_);
    partition->setLen.resize(g->numTot_);
    partition->setCount.resize(g->numTot_);
    partition->setFront.resize(g->numTot_);
    // Allocate refinement stuff
    refSet.resize(g->numTot_);
    splitSet.resize(g->numTot_);
    newSetSize.resize(g->numTot_ + 2);
    sList.assign(g->numTot_, -1);
    edgFreq.assign(g->numTot_, 0);
    vEdgColor.assign(g->nnz_, 0);
    vEdg.assign(g->nnz_, 0);
    wghtFreq.assign(g->numWeights_, 0);
    wghtStart.assign(g->numWeights_ + 1, 0);
    wghtEdg.assign(g->nnz_, 0);
    wghtOffset.assign(g->numWeights_, 0);
    scount.assign(g->numTot_, 0);
    count.assign(g->numTot_ + 1, 0);
    conncnts.assign(g->numTot_, 0);
    indMark.assign(g->numTot_, false);
    nInd.assign(g->numTot_, 0);
    sInd.assign(g->numTot_, 0);
    nextNon.assign(g->numTot_ + 1, 0);
    prevNon.assign(g->numTot_ + 1, 0);
    // Additions hopefully can further improve
    cDeg.assign(g->numTot_, 0);
    cDegFreq.assign(maxDeg + 1, 0);
    cDegIndex.assign(maxDeg + 1, 0);
    indexCDeg.assign(maxDeg + 1, 0);
    splitOffset.assign(maxDeg + 1, 0);
    newFront.assign(g->numTot_, 0);
    adjCell.assign(g->numTot_, 0);
    nodeAdj.assign(g->numTot_, 0);
    cellAdj.assign(g->numTot_, 0);
    refSize.assign(g->numTot_, 0);
    size.assign(g->numTot_, 0);
    // Fill in partition
    for (i = 0; i < g->numTot_; ++i){
        partition->setCount[g->colors[i]]++;
        if (max < g->colors[i]) max = g->colors[i]; 
    }
    partition->setSize = partition->setCount;
    partition->setLen[0] = partition->setCount[0] - 1;
    size[0] = partition->setLen[0] + 1;
    for (i = 0; i < max; ++i){
        partition->setLen[partition->setCount[i]] = partition->setCount[i + 1] - 1;
        partition->setCount[i + 1] += partition->setCount[i];
        size[partition->setCount[i]] = partition->setLen[partition->setCount[i]] + 1; 
    }
    for (i = 0; i < g->numTot_; ++i){
        setLabel(--partition->setCount[g->colors[i]], i);
        partition->set[i] = g->colors[i];
    }
    for (i = 0; i < g->numTot_; i += partition->setLen[i] + 1){
        addInduce(i);
        fixFronts(i, i);
    }
    // Prepare the nextNon and prevNon lists
    for (i = 0, j = -1; i < g->numTot_; i += partition->setLen[i] + 1){
        if (!partition->setLen[i]) continue;
        prevNon[i] = j;
        nextNon[j + 1] = i;
        j = i;
    }
    prevNon[g->numTot_] = j;
    nextNon[j] = g->numTot_;
    label.assign(partition->label.begin(), partition->label.end());
    junk.assign(g->numTot_, 0);
    front.assign(partition->front.begin(), partition->front.end());
    len.assign(partition->setLen.begin(), partition->setLen.end());
    refine();
}

void HighsOCEquitablePartition::lp2Graph(){
    // Ints and arrays for converting 
    int i, j, w, tempk, tempj, nColor = 0, nTot, nnz = originalLp->Avalue_.size(), 
    numCol = originalLp->numCol_, numRow = originalLp->numRow_;
    int *aout, *ain, *eout, *ein, *wout, *win, *colors, *edgColors;
    nTot = numCol + numRow;
    nnz = originalLp->Avalue_.size();
    aout = (int *)calloc( (nTot+1), sizeof(int) );
    eout = (int *)calloc( 2 * nnz, sizeof(int) );
    wout = (int *)calloc( 2 * nnz, sizeof(int) ); /* weight data */
    colors = (int *)calloc( nTot, sizeof(int) );
    edgColors = (int *)calloc( 2 *nnz, sizeof(int) );
    // Maps used to create {P}^0
    map<double, int> wghtColors;
    map<tuple<double, double>, int> rowColors;
	map<tuple<double, double, double>, int > colColors;
    g = (struct OCgraph*)calloc(1, sizeof(struct OCgraph));
    g->nnz_ = nnz;
    g->numCol_ = numCol;
    g->numRow_ = numRow;
    g->numTot_ = numCol + numRow;
    // Set g pointers correctly
    g->adj = aout;
    g->edg = eout;
    g->wght = wout;
    g->colors = colors;
    g->edgColors = edgColors;
    // Set in = out
    ain = aout;
	ein = eout;
	win = wout;
    // Fill in {P}^0 not necessarily equitable
    pair<map<tuple<double, double, double>, int>::iterator, bool> retCol;
    pair<map<tuple<double, double>, int>::iterator, bool> retRow;
    pair<map<double, int>::iterator, bool> retWght;
    // Partition columns based on lb, ub, c^t
    for (i = 0; i < numCol; ++i){
        retCol = colColors.insert(pair<tuple<double, double, double>, int>(
			make_tuple(originalLp->colLower_[i], originalLp->colUpper_[i], 
            originalLp->colCost_[i]), nColor));
        colors[i] = colColors.find(make_tuple(originalLp->colLower_[i], originalLp->colUpper_[i],
        originalLp->colCost_[i]))->second;
        if (retCol.second) ++nColor;
    }
    // Parititon rows based on lb, ub
    for (i = 0; i < numRow; ++i){
        retRow = rowColors.insert(pair<tuple<double, double>, int>(
			make_tuple(originalLp->rowLower_[i], originalLp->rowUpper_[i]), 
            nColor));
        colors[i + numCol] = rowColors.find(make_tuple(originalLp->rowLower_[i], originalLp->rowUpper_[i]))->second;
        if (retRow.second) ++nColor;
    }
    // Color matrix coefficients for easier computation
    nColor = 0;
    for (i = 0; i < nnz; ++i){
        retWght = wghtColors.insert(pair<double,int>(originalLp->Avalue_[i],
            nColor));
        if (retWght.second) ++nColor;
    }
    // Set graph num diffferent weights
    g->numWeights_ = wghtColors.size();
    // Fill in adj matrix
    for (j = 0; j < numCol; ++j){
        for (i = originalLp->Astart_[j]; i < originalLp->Astart_[j + 1]; ++i){
            ++ain[originalLp->Aindex_[i] + numCol]; ++aout[j];
        }
    }
    // Fix adj matrix
    fixAdjMiddle(nTot, aout);
    // Insert the actual wght data and connections
    for (j = 0; j < numCol; ++j){
        for (i = originalLp->Astart_[j]; i < originalLp->Astart_[j + 1]; ++i){
            tempk = ain[originalLp->Aindex_[i] + numCol]++; tempj = aout[j]++;
            eout[tempj] = originalLp->Aindex_[i] + numCol; ein[tempk] = j;
            w = wghtColors.find(originalLp->Avalue_[i])->second;
            wout[tempj] = w; win[tempk] = w;
            edgColors[tempj] = w; edgColors[tempk] = w;
        }
    }
    // Fix the adjacency ends
    fixAdjEnds(nTot, 2 * nnz, aout);
    int sum;
    maxDeg = 0;
    for (j = 0; j < nTot; ++j){
        sum = 0;
        for (i = g->adj[j]; i < g->adj[j + 1]; ++i)
            sum += g->wght[i] + 1;
        if (sum > maxDeg) maxDeg = sum;
    }
}

void HighsOCEquitablePartition::fixAdjMiddle(int n, int* adj){
    int val, sum, i;
    val = adj[0]; sum = 0; adj[0] = 0;
    for (i = 1; i < n; ++i){
        sum += val;
        val = adj[i];
        adj[i] = sum;
    }
}

void HighsOCEquitablePartition::fixAdjEnds(int n, int e, int* adj){
    int i;
    for (i = n - 1; i > 0; --i)
        adj[i] = adj[i - 1];
    adj[0] = 0;
    adj[n] = e;
}

void HighsOCEquitablePartition::fixFronts(int cf, int ff){
    /* cf, ff are current front and fixed front resp. */
    int i, end = cf + partition->setLen[cf];
    for (i = ff; i <= end; ++i){
        partition->front[partition->label[i]] = cf;
    }
}

void HighsOCEquitablePartition::setLabel(int index, int value){
    partition->label[index] = value;
    partition->unlabel[value] = index;
}

void HighsOCEquitablePartition::dataCount(int edg){
    int cf = partition->front[edg];
    if (partition->setLen[cf] && !scount[edg]++) moveTo(edg);
}

void HighsOCEquitablePartition::dataMark(int edg){
    int cf = partition->front[edg];
    if (partition->setLen[cf]) moveTo(edg);
}

void HighsOCEquitablePartition::moveTo(int edg){
    int cf = partition->front[edg];
    int cb = cf + partition->setLen[cf];
    int offset = conncnts[cf]++;
    swapLabels(cb - offset, partition->unlabel[edg]);
    if (!offset) sList[sListSize++] = cf;
}

void HighsOCEquitablePartition::swapLabels(int a, int b){
    int tmp = partition->label[a];
    setLabel(a, partition->label[b]);
    setLabel(b, tmp);
}

void HighsOCEquitablePartition::addInduce(int cf){
    if (!partition->setLen[cf])
        sInd[sIndSize++] = cf;
    else
        nInd[nIndSize++] = cf;
    indMark[cf] = true;
    S.push(cf);
}

bool HighsOCEquitablePartition::discrete(){
    return nSplits == g->numTot_;
}

void HighsOCEquitablePartition::clear(){
    int i = 0; 
    for (i = 0; i < nIndSize; ++i)
        indMark[nInd[i]] = false;
    for (i = 0; i < sIndSize; ++i)
        indMark[sInd[i]] = false;
    nIndSize = sIndSize = 0;
}

void HighsOCEquitablePartition::markCell(int cf){
    adjCell[adjCellCnt++] = cf;
}

void HighsOCEquitablePartition::splitCell(int cf){
    int i, j = label[cf], k, edg, cb = cf + len[cf], 
    nf = cf, cdeg, numCDeg = 0, lab, idx;
    /* if size[cf] != refSize[cf], numCDeg += 1 */
    if (refSize[cf] != size[cf]){
        cDegIndex[0] = numCDeg;
        indexCDeg[numCDeg++] = 0;
        cDegFreq[0] = size[cf] - refSize[cf];
    }
    for (k = cf; k < cb + 1; ++k){
        edg = label[k];
        if(cDeg[edg] && !cDegFreq[cDeg[edg]]++){
            cDegIndex[cDeg[edg]] = numCDeg;
            indexCDeg[numCDeg++] = cDeg[edg];
        }
    }
    newFront[0] = nf;
    for (k = 1; k < numCDeg; ++k){
        cdeg = indexCDeg[k - 1];
        splitOffset[k] = cDegFreq[cdeg];
        nf += splitOffset[k];
        newFront[k] = nf;
    }
    for (k = cf; k < cb + 1; ++k){
        lab = label[k];
        cdeg = cDeg[lab];
        idx = cDegIndex[cdeg];
        nf = cf + splitOffset[idx]++;
        junk[nf] = lab;
    }
    for (k = cf; k < cb + 1; ++k){
        cdeg = cDeg[junk[k]];
        idx = cDegIndex[cdeg];
        label[k] = junk[k];
        front[k] = newFront[idx];
    }
    // for (i = 0; i < numCDeg; ++i){
    //     newFront[i] = nf;
    //     nf += cDegFreq[]
    // }
    for (k = 0; k < numCDeg; ++k){
        cdeg = indexCDeg[k];
        splitOffset[k] = 0;
        cDegIndex[k] = 0;
        cDegFreq[cdeg] = 0;
    }
    
}