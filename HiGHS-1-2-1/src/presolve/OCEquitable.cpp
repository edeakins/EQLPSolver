#include "OCEquitable.h"
using namespace std;

OCPartition HighsOCEquitablePartition::getPartition(){
    return *partition;
}

void HighsOCEquitablePartition::runToDiscrete(){
    while (!discrete()) isolate();
}

bool HighsOCEquitablePartition::isolate(){
    // int i, targ;
    // for (i = 0; i < g->numTot_; i += len[i] + 1){
    //     if (len[i]){
    //         targ = i;
    //         break;
    //     }
    // }
    partition->level++;
    int targ = nextNon[-1];
    int min = targ;
    // int length = len[targ];
    int back = targ + partition->len[targ];
    int lab = partition->label[back];
    // int nsf = back;
    // // Additions, hopefully better
    // front[back] = nsf; 
    // index[lab] = front[nsf]; 
    // len[targ]--; len[nsf] = 0;
    // size[targ]--; size[nsf] = 1; ++nSplits;
    // addInduce(nsf);
    swapLabels(min, back);
    split(targ, back);
    refine();
    // for (int i = 0; i < g->numTot_; i += partition->len[i] + 1)
    //     sort(partition->label.begin() + i, partition->label.begin() + i + partition->len[i] + 1);
    if (discrete()) return true;
    else return false;
}

bool HighsOCEquitablePartition::refine(){
    int f;
    while (true){
        if (discrete()){ 
            clear();
            return true;
        }
        if (sIndSize){
            f = sInd[--sIndSize];
            indMark[f] = false;
            if (!refineSingleCell(f)) break;
        }
        else if (nIndSize){ 
            f = nInd[--nIndSize];
            indMark[f] = false;
            if (!refineNonsingleCell(f)) break;
        }
        else{
            return true;
        }
    }
    clear();
    return false;
}

// Refining cell was nonsingleton
bool HighsOCEquitablePartition::refineNonsingleCell(int sf){
    int i, j, k, edg, wght, sb, cf, ssize, e, w, 
    wcount = 0, ecount = 0, numAdj = 0, idx;
    bool ret = true;
    sb = sf + partition->len[sf];
    /* Double check for nonsingles which became singles later */
    if (sf == sb) {
        return refineSingleCell(sf);
    }
    // sb = sf + len[sf];
    // maxDeg = 0; adjCellCnt = 0;
    // Compute degree sums for all vertices adjacent to refining cell
    for (i = sf; i < sb + 1; ++i){
        j = partition->label[i];
        // j = label[i];
        for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
            edg = g->edg[k]; wght = g->wght[k]; //cf = partition->front[edg];
            // edg = g->edg[k]; wght = g->wght[k]; cf = front[index[edg]];
            if (!wghtFreq[wght]++) wcountIdx[wcount++] = wght;
            // vEdgColor[ecount] = wght;
            // vEdg[ecount++] = edg;
            wghtCnt[wght].push_back(edg);
            // // Additions to move away from saucy
            // markCell(cf, edg, wght);
        }
    }

    // // Additions to move away from saucy
    // for (k = 0; k < adjCellCnt; ++k){
    //     splitCell(adjCell[k]);
    //     refSize[adjCell[k]] = 0;
    //     cellAdj[adjCell[k]] = 0;
    //     adjCell[k] = 0;
    // }
    // for (i = sf; i < sb + 1; ++i){
    //     // j = partition->label[i];
    //     j = label[i];
    //     for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
    //         edg = g->edg[k];
    //         if (cDeg[edg]) cDeg[edg] = 0;
    //         if (nodeAdj[edg]) nodeAdj[edg] = 0;
    //     }
    // }
    // // Populate sparse start array for edge weights
    // for (i = 0; i < wcount; ++i){
    //     wght = wcountIdx[i];
    //     wghtStart[i + 1] += (wghtStart[i] + wghtFreq[wght]); // here fix this
    //     wghtOffset[wght] = wghtStart[i];
    // }

    // // Populate sparse weight index array.  
    // // Don't need sparse val array as the indices are stored
    // // according to their weight.
    // for (i = 0; i < ecount; ++i){
    //     w = vEdgColor[i];
    //     e = vEdg[i];
    //     idx = wghtOffset[w]++;
    //     wghtEdg[idx] = e;
    // }

    for (i = 0; i < wcount; ++i){
        wght = wcountIdx[i];
        wcountIdx[i] = 0;
        for (j = 0; j < wghtCnt[wght].size(); ++j){
            dataCount(wghtCnt[wght][j]);
        }
        ret = ret && refineNonSingles();
        for (j = 0; j < wghtCnt[wght].size(); ++j){
            scount[wghtCnt[wght][j]] = 0;
        }
        // for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
        //     scount[wghtEdg[j]] = 0;
        //     wghtEdg[j] = 0;
        //     vEdgColor[j] = 0;
        //     vEdg[j] = 0;
        // }
        // wghtStart[i] = 0;
        // wghtOffset[wght] = 0;
        wghtFreq[wght] = 0;
        wghtCnt[wght].clear();
    }
    // wghtStart[wcount] = 0;
    return ret;
}

// Refining cell is a singleton
bool HighsOCEquitablePartition::refineSingleCell(int sf){
    int i, j = partition->label[sf], k, edg, wght, sb, cf, 
    ssize, e, w, wcount = 0, ecount = 0, numAdj = 0, idx, cdeg;
    bool ret = true;
    // maxDeg = 0; adjCellCnt = 0;
    for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
        // edg = g->edg[k]; wght = g->wght[k]; cf = partition->front[edg];
        edg = g->edg[k]; wght = g->wght[k];// cf = front[index[edg]];
        if (!wghtFreq[wght]++) wcountIdx[wcount++] = wght;
        // vEdgColor[ecount] = wght;
        // vEdg[ecount++] = edg;
        wghtCnt[wght].push_back(edg);
        // // Additions to move away from saucy
        // markCell(cf, edg, wght);
    }

    // // Additions to move away from saucy
    // for (k = 0; k < adjCellCnt; ++k){
    //     splitCell(adjCell[k]);
    //     refSize[adjCell[k]] = 0;
    //     cellAdj[adjCell[k]] = 0;
    //     adjCell[k] = 0;
    // }
    // for (k = g->adj[j]; k < g->adj[j + 1]; ++k){
    //     edg = g->edg[k];
    //     if (cDeg[edg]) cDeg[edg] = 0;
    //     if (nodeAdj[edg]) nodeAdj[edg] = 0;
    // }
    // // Populate sparse start array for edge weights
    // for (i = 0; i < wcount; ++i){
    //     wght = wcountIdx[i];
    //     wghtStart[i + 1] += (wghtStart[i] + wghtFreq[wght]);
    //     wghtOffset[wght] = wghtStart[i];
    // }

    // // Populate sparse weight index array.  
    // // Don't need sparse val array as the indices are stored
    // // according to their weight.
    // for (i = 0; i < ecount; ++i){
    //     w = vEdgColor[i];
    //     e = vEdg[i];
    //     idx = wghtOffset[w]++;
    //     wghtEdg[idx] = e;
    // }

    for (i = 0; i < wcount; ++i){
        wght = wcountIdx[i];
        wcountIdx[i] = 0;
        for (j = 0; j < wghtCnt[wght].size(); ++j){
            dataMark(wghtCnt[wght][j]);
        }
        ret = ret && refineNonSingles();
        for (j = 0; j < wghtCnt[wght].size(); ++j){
            scount[wghtCnt[wght][j]] = 0;
        }
        // for (j = wghtStart[i]; j < wghtStart[i + 1]; ++j){
        //     scount[wghtEdg[j]] = 0;
        //     wghtEdg[j] = 0;
        //     vEdgColor[j] = 0;
        //     vEdg[j] = 0;
        // }
        // wghtStart[i] = 0;
        // wghtOffset[wght] = 0;
        wghtFreq[wght] = 0;
        wghtCnt[wght].clear();
    }
    // wghtStart[wcount] = 0;
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
    cb = cf + partition->len[cf];
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
    int zcnt = partition->len[cf] + 1 - conncnts[cf];
    return possiblySplit(cf, cf + zcnt);
}

// Does the acutal splitting of a cell once new cell sizes are known
bool HighsOCEquitablePartition::split(int cf, int ff){
    // Track number of splits
    ++partition->nsplits;
    cf < g->numCol_ ? ++partition->ncsplits : ++partition->nrsplits;
    // Do splitting
    int cb, fb;
    fb = ff - 1;
    cb = cf + partition->len[cf];
    partition->len[cf] = fb - cf;
    partition->len[ff] = cb - ff;
    // Fix the fronts
    fixFronts(ff, ff);
    // Add inducers
    if (indMark[cf] || partition->len[ff] < partition->len[cf])
        addInduce(ff);
    else   
        addInduce(cf);
    // Keep track of non singletons for new targets for inducing another level of refinement
    if (partition->len[ff]){
        prevNon[nextNon[cf]] = ff;
        nextNon[ff] = nextNon[cf];
        prevNon[ff] = cf;
        nextNon[cf] = ff;
    }
    if (!partition->len[cf]){
        nextNon[prevNon[cf]] = nextNon[cf];
        prevNon[nextNon[cf]] = prevNon[cf];
    }
    return true;
}

// Possibly does actual splitting when there is a 0 wght frequency
bool HighsOCEquitablePartition::possiblySplit(int cf, int ff){
    return cf == ff ? true : split(cf, ff);
}

bool HighsOCEquitablePartition::allocatePartition(){
    // Convert LP to sparse graph rep
    lp2Graph();
    int i, j, max = 0, maxSetSize = 0, idx = 0;
    // maxSetSize = g->numCol_ > g->numRow_ ? g->numCol_ : g->numRow_;
    // Allocate partition storage
    partition = (OCPartition*)calloc(1, sizeof(OCPartition));
    partition->target = -1;
    partition->level = 0;
    partition->nsplits = 0;
    partition->ncsplits = 0;
    partition->nrsplits = 0;
    partition->front.resize(g->numTot_);
    partition->label.resize(g->numTot_);
    partition->unlabel.resize(g->numTot_);
    partition->len.resize(g->numTot_);
    // partition->parent.resize(g->numTot_);
    // partition->set.resize(g->numTot_);
    // partition->setSize.resize(g->numTot_);
    // partition->setFront.resize(g->numTot_);
    // Allocate refinement stuff
    refSet.resize(g->numTot_);
    splitSet.resize(g->numTot_);
    newSetSize.resize(g->numTot_ + 2);
    sList.assign(g->numTot_, 0);
    // edgFreq.assign(g->numTot_, 0);
    // vEdgColor.assign(g->nnz_, 0);
    // vEdg.assign(g->nnz_, 0);
    wghtFreq.assign(g->numWeights_, 0);
    wcountIdx.assign(g->numWeights_, 0);
    // wghtStart.assign(g->numWeights_ + 1, 0);
    // wghtEdg.assign(g->nnz_, 0);
    // wghtOffset.assign(g->numWeights_, 0);
    scount.assign(g->numTot_, 0);
    count.assign(g->numTot_ + 1, 0);
    conncnts.assign(g->numTot_, 0);
    indMark.assign(g->numTot_, false);
    nInd.assign(g->numTot_, 0);
    sInd.assign(g->numTot_, 0);
    nextNon = allocInts(g->numTot_ + 1) + 1;
    prevNon = allocInts(g->numTot_ + 1);
    // // Additions hopefully can further improve
    // cDeg.assign(g->numTot_, 0);
    // cDegFreq.assign(maxDeg + 1, 0);
    // cDegIndex.assign(maxDeg + 1, 0);
    // indexCDeg.assign(maxDeg + 1, 0);
    // splitOffset.assign(maxDeg + 1, 0);
    // newFront.assign(g->numTot_, 0);
    // adjCell.assign(g->numTot_, 0);
    // nodeAdj.assign(g->numTot_, 0);
    // cellAdj.assign(g->numTot_, 0);
    // refSize.assign(g->numTot_, 0);
    // label.assign(g->numTot_, 0);
    // front.assign(g->numTot_, 0);
    // index.assign(g->numTot_, 0);
    // len.assign(g->numTot_, 0);
    // size.assign(g->numTot_, 0);
    // set.assign(g->numTot_, 0);
    // junk.assign(g->numTot_, 0);
    // Fill in partition
    for (i = 0; i < g->numTot_; ++i){
        scount[g->colors[i]]++;
        // count[g->colors[i]]++;
        if (max < g->colors[i]) max = g->colors[i]; 
    }
    // partition->setSize = scount;
    partition->len[0] = scount[0] - 1;
    // len[0] = count[0] - 1;
    // size[0] = len[0] + 1;
    for (i = 0; i < max; ++i){
        partition->len[scount[i]] = scount[i + 1] - 1;
        scount[i + 1] += scount[i];
        // partition->setSize[scount[i]] = partition->len[scount[i]] + 1;
        // len[count[i]] = count[i + 1] - 1;
        // count[i + 1] += count[i];
        // size[count[i]] = len[count[i]] + 1; 
    }
    partition->nsplits = max + 1;
    for (i = 0; i < g->numTot_; ++i){
        setLabel(--scount[g->colors[i]], i);
        // setLabel(--count[g->colors[i]], i);
        // set[i] = g->colors[i];
    }
    for (i = 0; i < g->numTot_; i += partition->len[i] + 1){
        addInduce(i);
        fixFronts(i, i);
        partition->front[i] < g->numCol_ ? partition->ncsplits++ : partition->nrsplits++;
    }
    // Prepare the nextNon and prevNon lists
    for (i = 0, j = -1; i < g->numTot_; i += partition->len[i] + 1){
        if (!partition->len[i]) continue;
        prevNon[i] = j;
        nextNon[j] = i;
        j = i;
    }
    prevNon[g->numTot_] = j;
    nextNon[j] = g->numTot_;
    /* Clear out scount */
    for (i = 0; i <= max; ++i)
        scount[i] = 0;
    // label.assign(partition->label.begin(), partition->label.end());
    // front.assign(partition->front.begin(), partition->front.end());
    // len.assign(partition->len.begin(), partition->len.end());
    refine();
    if (discrete()) return true;
    return false;
}

bool HighsOCEquitablePartition::allocatePartition(HighsLp* lp){
    // Convert LP to sparse graph rep
    originalLp = lp;
    lp2Graph();
    int i, j, max = 0, maxSetSize = 0, idx = 0;
    // maxSetSize = g->numCol_ > g->numRow_ ? g->numCol_ : g->numRow_;
    // Allocate partition storage
     partition = (OCPartition*)calloc(1, sizeof(OCPartition));
    partition->target = -1;
    partition->level = 0;
    partition->nsplits = 0;
    partition->ncsplits = 0;
    partition->nrsplits = 0;
    partition->front.resize(g->numTot_);
    partition->label.resize(g->numTot_);
    partition->unlabel.resize(g->numTot_);
    partition->len.resize(g->numTot_);
    // partition->parent.resize(g->numTot_);
    // partition->set.resize(g->numTot_);
    // partition->setSize.resize(g->numTot_);
    // partition->setFront.resize(g->numTot_);
    // Allocate refinement stuff
    refSet.resize(g->numTot_);
    splitSet.resize(g->numTot_);
    newSetSize.resize(g->numTot_ + 2);
    sList.assign(g->numTot_, 0);
    // edgFreq.assign(g->numTot_, 0);
    // vEdgColor.assign(g->nnz_, 0);
    // vEdg.assign(g->nnz_, 0);
    wghtFreq.assign(g->numWeights_, 0);
    wcountIdx.assign(g->numWeights_, 0);
    // wghtStart.assign(g->numWeights_ + 1, 0);
    // wghtEdg.assign(g->nnz_, 0);
    // wghtOffset.assign(g->numWeights_, 0);
    scount.assign(g->numTot_, 0);
    count.assign(g->numTot_ + 1, 0);
    conncnts.assign(g->numTot_, 0);
    indMark.assign(g->numTot_, false);
    nInd.assign(g->numTot_, 0);
    sInd.assign(g->numTot_, 0);
    nextNon = allocInts(g->numTot_ + 1) + 1;
    prevNon = allocInts(g->numTot_ + 1);
    // // Additions hopefully can further improve
    // cDeg.assign(g->numTot_, 0);
    // cDegFreq.assign(maxDeg + 1, 0);
    // cDegIndex.assign(maxDeg + 1, 0);
    // indexCDeg.assign(maxDeg + 1, 0);
    // splitOffset.assign(maxDeg + 1, 0);
    // newFront.assign(g->numTot_, 0);
    // adjCell.assign(g->numTot_, 0);
    // nodeAdj.assign(g->numTot_, 0);
    // cellAdj.assign(g->numTot_, 0);
    // refSize.assign(g->numTot_, 0);
    // label.assign(g->numTot_, 0);
    // front.assign(g->numTot_, 0);
    // index.assign(g->numTot_, 0);
    // len.assign(g->numTot_, 0);
    // size.assign(g->numTot_, 0);
    // set.assign(g->numTot_, 0);
    // junk.assign(g->numTot_, 0);
    // Fill in partition
    for (i = 0; i < g->numTot_; ++i){
        scount[g->colors[i]]++;
        // count[g->colors[i]]++;
        if (max < g->colors[i]) max = g->colors[i]; 
    }
    // partition->setSize = scount;
    partition->len[0] = scount[0] - 1;
    // len[0] = count[0] - 1;
    // size[0] = len[0] + 1;
    for (i = 0; i < max; ++i){
        partition->len[scount[i]] = scount[i + 1] - 1;
        scount[i + 1] += scount[i];
        // partition->setSize[scount[i]] = partition->len[scount[i]] + 1;
        // len[count[i]] = count[i + 1] - 1;
        // count[i + 1] += count[i];
        // size[count[i]] = len[count[i]] + 1; 
    }
    partition->nsplits = max + 1;
    for (i = 0; i < g->numTot_; ++i){
        setLabel(--scount[g->colors[i]], i);
        // setLabel(--count[g->colors[i]], i);
        // set[i] = g->colors[i];
    }
    for (i = 0; i < g->numTot_; i += partition->len[i] + 1){
        addInduce(i);
        fixFronts(i, i);
        partition->front[i] < g->numCol_ ? partition->ncsplits++ : partition->nrsplits++;
    }
    // Prepare the nextNon and prevNon lists
    for (i = 0, j = -1; i < g->numTot_; i += partition->len[i] + 1){
        if (!partition->len[i]) continue;
        prevNon[i] = j;
        nextNon[j] = i;
        j = i;
    }
    prevNon[g->numTot_] = j;
    nextNon[j] = g->numTot_;
    for (i = 0; i <= max; ++i)
        scount[i] = 0;
    // label.assign(partition->label.begin(), partition->label.end());
    // front.assign(partition->front.begin(), partition->front.end());
    // len.assign(partition->len.begin(), partition->len.end());
    refine();
    // for (int i = 0; i < g->numTot_; i += partition->len[i] + 1)
    //     sort(partition->label.begin() + i, partition->label.begin() + i + partition->len[i] + 1);
    if (discrete()) return true;
    return false;
}

void HighsOCEquitablePartition::lp2Graph(){
    // Ints and arrays for converting 
    int i, j, w, tempk, tempj, nColor = 0, nTot, nnz = originalLp->a_matrix_.value_.size(), 
    numCol = originalLp->num_col_, numRow = originalLp->num_row_;
    int *aout, *ain, *eout, *ein, *wout, *win, *colors, *edgColors;
    nTot = numCol + numRow;
    nnz = originalLp->a_matrix_.value_.size();
    aout = (int *)calloc( (nTot+1), sizeof(int) );
    eout = (int *)calloc( 2 * nnz, sizeof(int) );
    wout = (int *)calloc( 2 * nnz, sizeof(int) ); /* weight data */
    colors = (int *)calloc( nTot, sizeof(int) );
    edgColors = (int *)calloc( 2 *nnz, sizeof(int) );
    // Maps used to create {P}^0
    map<double, int> wghtColors;
    map<tuple<double, double>, int> rowColors;
	map<tuple<double, double, double>, int > colColors;
    g = (OCGraph*)calloc(1, sizeof(OCGraph));
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
			make_tuple(originalLp->col_lower_[i], originalLp->col_upper_[i], 
            originalLp->col_cost_[i]), nColor));
        colors[i] = colColors.find(make_tuple(originalLp->col_lower_[i], originalLp->col_upper_[i],
        originalLp->col_cost_[i]))->second;
        if (retCol.second) ++nColor;
    }
    // Parititon rows based on lb, ub
    for (i = 0; i < numRow; ++i){
        retRow = rowColors.insert(pair<tuple<double, double>, int>(
			make_tuple(originalLp->row_lower_[i], originalLp->row_upper_[i]), 
            nColor));
        colors[i + numCol] = rowColors.find(make_tuple(originalLp->row_lower_[i], originalLp->row_upper_[i]))->second;
        if (retRow.second) ++nColor;
    }
    // Color matrix coefficients for easier computation
    nColor = 0;
    for (i = 0; i < nnz; ++i){
        retWght = wghtColors.insert(pair<double,int>(originalLp->a_matrix_.value_[i],
            nColor));
        if (retWght.second) ++nColor;
    }
    // Set graph num diffferent weights
    g->numWeights_ = wghtColors.size();
    // Fill in adj matrix
    for (j = 0; j < numCol; ++j){
        for (i = originalLp->a_matrix_.start_[j]; i < originalLp->a_matrix_.start_[j + 1]; ++i){
            ++ain[originalLp->a_matrix_.index_[i] + numCol]; ++aout[j];
        }
    }
    // Fix adj matrix
    fixAdjMiddle(nTot, aout);
    // Insert the actual wght data and connections
    for (j = 0; j < numCol; ++j){
        for (i = originalLp->a_matrix_.start_[j]; i < originalLp->a_matrix_.start_[j + 1]; ++i){
            tempk = ain[originalLp->a_matrix_.index_[i] + numCol]++; tempj = aout[j]++;
            eout[tempj] = originalLp->a_matrix_.index_[i] + numCol; ein[tempk] = j;
            w = wghtColors.find(originalLp->a_matrix_.value_[i])->second;
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
    int i, end = cf + partition->len[cf];
    // int i, end = cf + len[cf];
    for (i = ff; i <= end; ++i){
        partition->front[partition->label[i]] = cf;
        // front[label[i]] = cf;
    }
}

void HighsOCEquitablePartition::setLabel(int idx, int val){
    partition->label[idx] = val;
    partition->unlabel[val] = idx;
    // label[idx] = val;
    // index[val] = idx;
}

void HighsOCEquitablePartition::dataCount(int edg){
    int cf = partition->front[edg];
    if (partition->len[cf] && !scount[edg]++) moveTo(edg);
}

void HighsOCEquitablePartition::dataMark(int edg){
    int cf = partition->front[edg];
    if (partition->len[cf]) moveTo(edg);
}

void HighsOCEquitablePartition::moveTo(int edg){
    int cf = partition->front[edg];
    int cb = cf + partition->len[cf];
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
    if (!partition->len[cf])
        sInd[sIndSize++] = cf;
    else
        nInd[nIndSize++] = cf;
    indMark[cf] = true;
    // S.push(cf);
}

bool HighsOCEquitablePartition::discrete(){
    // if (partition->ncsplits == numColSets &&
    //     partition->nrsplits == numRowSets)
    //     return true;
    // else{
    //     numColSets = partition->ncsplits;
    //     numRowSets = partition->nrsplits;
    // }
    // return partition->nsplits == g->numTot_ ||
    //        partition->ncsplits == g->numCol_ ||
    //        partition->nrsplits == g->numRow_;
    return partition->nsplits == g->numTot_;
}

void HighsOCEquitablePartition::clear(){
    int i = 0; 
    for (i = 0; i < nIndSize; ++i)
        indMark[nInd[i]] = false;
    for (i = 0; i < sIndSize; ++i)
        indMark[sInd[i]] = false;
    nIndSize = sIndSize = 0;
}

void HighsOCEquitablePartition::markCell(int cf, int edg, int wght){
    if ((cDeg[edg] += wght + 1) > maxDeg) maxDeg = cDeg[edg];
    if (len[cf]  && !cellAdj[cf]){
        cellAdj[cf]++;
        adjCell[adjCellCnt++] = cf;
    }
    if (len[cf] && !nodeAdj[edg]){
        nodeAdj[edg]++;
        refSize[cf]++;
    }
}

int * HighsOCEquitablePartition::allocInts(int n){
    int *p = (int *)malloc(n * sizeof(int));
    return p;
}

void HighsOCEquitablePartition::splitCell(int cf){
    int i, j = label[cf], k, edg, cb = cf + len[cf], 
    nf = cf, cdeg, numCDeg = 0, lab, idx, maxSSize, maxSizeFront;
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
    if (numCDeg == 1){
        for (k = 0; k < numCDeg; ++k){
            cdeg = indexCDeg[k];
            cDegIndex[k] = 0;
            cDegFreq[cdeg] = 0;
        }
        return;
    }
    splitOffset[0] = nf;
    newFront[0] = nf;
    len[nf] = cDegFreq[indexCDeg[0]] - 1;
    size[nf] = cDegFreq[indexCDeg[0]];
    maxSSize = size[nf]; maxSizeFront = nf;
    for (k = 1; k < numCDeg; ++k){
        cdeg = indexCDeg[k - 1];
        splitOffset[k] = newFront[k - 1] + cDegFreq[cdeg];
        nf = splitOffset[k];
        newFront[k] = nf;
        cdeg = indexCDeg[k];
        len[nf] = cDegFreq[cdeg] - 1;
        if ((size[nf] = cDegFreq[cdeg]) > maxSSize){
            maxSSize = size[nf];
            maxSizeFront = nf;
        } 
        ++nSplits;
    }
    for (k = cf; k < cb + 1; ++k){
        lab = label[k];
        cdeg = cDeg[lab];
        idx = cDegIndex[cdeg];
        nf = splitOffset[idx]++;
        junk[nf] = lab;
        index[lab] = nf;
    }
    for (k = cf; k < cb + 1; ++k){
        cdeg = cDeg[junk[k]];
        idx = cDegIndex[cdeg];
        label[k] = junk[k];
        front[k] = newFront[idx];
        if (!(front[k] == maxSizeFront) && !indMark[front[k]]) addInduce(front[k]);
    }
    for (k = 0; k < numCDeg; ++k){
        cdeg = indexCDeg[k];
        splitOffset[k] = 0;
        cDegIndex[k] = 0;
        cDegFreq[cdeg] = 0;
    }
}