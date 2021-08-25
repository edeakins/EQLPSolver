#ifndef PRESOLVE_OCEQUITABLE_H
#define PRESOLVE_OCEQUITABLE_H

// Includes
#include <map>
#include <string>
#include <stack>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
// #include <limits>
#include "HighsLp.h"

/* Aggregator class for LPs based off a equitable partition */
class HighsOCEquitablePartition{
public:
    /* Initialize and copy lp to equitable lp container */
    HighsOCEquitablePartition(){}
    HighsOCEquitablePartition(HighsLp* lp):
    originalLp(lp){
        allocatePartition();
    }
    /* //////////////////////////////////////////////////////////////
    SAUCY INSPIRED 
    ///////////////////////////////////////////////////////////// */
    /* If necessary, call run to discrete without intermediary code */
    void runToDiscrete();
    /* Calls isolation of nonsingleton class and new refinement round */
    void isolate();
    /* Refinement algorithm Berkholz 2017 */
    void refine();
    /* Refine cell */
    bool refineNonsingleCell(int sf);
    /* Refine cell based on sinlgeton */
    bool refineSinlgeCell(int sf);
    /* Calls splitting technique for all cells connected to nonsingle refining cell */
    bool refineNonSingles();
    /* Calls splitting technique for all cells connected to singleton refining cell */
    bool refineSingles();
    /* Split a refined cell */
    bool splitNonsingleCell(int cf);
    /* Split a refined cell when refining cell was singleton */
    bool splitSinlgeCell(int cf);
    /* Used to split the actual cells */
    bool split(int cf, int ff);
    /* Used to maybe split when there are 0 wght frequencies */
    bool possiblySplit(int cf, int ff);
    /* Allocate for all partition stuff */
    void allocatePartition();
    void allocatePartition(HighsLp* lp);
    /* Translate lp into bipartite graph */
    void lp2Graph();
    /* Fix adj middle positions */
    void fixAdjMiddle(int n, int* adj);
    /* Fix adj end positions */
    void fixAdjEnds(int n, int e, int* adj);
    /* Fix set fronts in partition storage */
    void fixFronts(int cf, int ff);
    /* Set the labels within sets of P */
    void setLabel(int index, int value);
    /* Mark cells for refinement when refining with a nonsingleton cell */
    void dataCount(int edge);
    /* Mark cells for refinement when refining with a singleton cell */
    void dataMark(int edge);
    /* Move elements to back of there cell */
    void moveTo(int edg);
    /* Swap labels around in label and unlabel array */
    void swapLabels(int a, int b);
    /* Add a cell front to induce list */
    void addInduce(int cf);
    /* Test whether partition is discrete */
    bool discrete();
    /* clear out induce info at the end of a ep */
    void clear();
    /* //////////////////////////////////////////////////////////////
    BERKHOLZ INSPIRED 
    ///////////////////////////////////////////////////////////// */
    /* Mark cells that are adjacent to refinement cell */
    void markCell(int cf, int edg, int wght);
    /* Split cells that were marked */
    void splitCell(int cf);


    /* HighsLp container for the original lp being passed in */
    HighsLp* originalLp;
    /* OCgraph storage for the lp for equitable calculations */
    OCgraph* g;
    /* Partition storage for the ep */
    OCpartition* partition;
    /* Refinement storage */
    int numColSets = 0;
    int numRowSets = 0;
    int sListSize= 0; // Number of sets in setList
    int nSplits = 0; // Number of splits that happen during refinement
    int nIndSize = 0; // Number of current nonsingleton inducers
    int sIndSize = 0; // Number of current singleton inducers
    std::stack<int> S; // Stack for refinement cells
    std::vector<int> refSet; // Contains set used for refinement
    std::vector<int> splitSet; // Contains the set to be split
    std::vector<int> newSetSize; // Contains the sizes of new sets from splitting
    std::vector<int> sList; // List of cells marked for refinement
    std::vector<int> edgFreq; 
    std::vector<int> vEdgColor; // Color of edg (u, v)
    std::vector<int> vEdg; // Keeps track of the current adj edges tor refining cell
    std::vector<int> wghtFreq; // Number of adjacencies to refinment cell with wght w
    std::vector<int> wghtOffset; // Keeps track of where next edg is places in wghtEdg and wght Index
    std::vector<int> wghtStart; // Sparse storage for refinement wghts
    std::vector<int> wghtEdg; // Sparse storage for refinement wghts
    std::vector<int> scount; // Number of connections to refining cell
     // Number of vertices with same adj count to ref cell
    std::vector<int> conncnts; // Connection counts for cell fronts
    std::vector<int> nInd; // Nonsingleton inducers of refinement
    std::vector<int> sInd; // Singleton inducers of refinement
    std::vector<bool> indMark; // Marks if a cell front is marked for refinement
    std::vector<int> nextNon; // Next non singleton cell from cf;
    std::vector<int> prevNon; // Prvious non singleton cell from cf;

    /*
     NOT CURRENTLY BEING USED
     */
    /* FOR .h */
    int maxDeg;
    int adjCellCnt;
    int numCDeg;
    // int maxCDeg;
    std::vector<int> cDeg;
    std::vector<int> cDegFreq;
    // std::vector<int> cDegStart;
    std::vector<int> cDegIndex;
    std::vector<int> indexCDeg;
    std::vector<int> splitOffset;
    std::vector<int> newFront;
    std::vector<int> adjCell;
    std::vector<int> cellAdj;
    std::vector<int> set;
    std::vector<int> label;
    std::vector<int> index;
    std::vector<int> junk;
    std::vector<int> front;
    std::vector<int> len;
    std::vector<int> size;
    std::vector<int> count;
    std::vector<int> refSize;
    std::vector<int> nodeAdj;
    // std::vector<bool> cDegB;
    // std::vector<int> cAdj;
    // std::vector<bool> cAdjB;
    // std::vector<int> maxCDeg;
    // std::vector<int> minCDeg;
    // std::vector<int> tempSet;
    // std::vector<int> tempSetSize;
    /* FOR .cpp */
    // cDeg.assign(g->numTot_, 0);
    // cDegStart.assign(g->numWeights_ + 1, 0);
    // cDegIndex.assign(g->nnz_, 0);
    // splitOffset.assign(maxCDeg, 0);
    // cDegB.assign(g->numTot_, false);
    // cDegFreq.assign(maxDeg, 0);
    // maxCDeg.assign(g->numTot_, 0);
    // minCDeg.assign(g->numTot_, kHighsIInf);
    // cAdj.assign(g->numTot_, -1);
    // cAdjB.assign(g->numTot_, false);
    // tempSet.assign(g->numTot_, -1);
    // tempSetSize.assign(g->numTot_, 0);
};

#endif