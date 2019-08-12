#ifndef HMODEL_H_
#define HMODEL_H_

#include "HMatrix.h"
#include "HFactor.h"
#include "HVector.h"
#include "HRandom.h"
#include "HTimer.h"

#include <string>
#include <vector>
#include <stack>
#include <set>
#include <list>
#include <map>
#include <tuple>
using namespace std;

const int LP_Status_Unset = -1;
const int LP_Status_Optimal = 0;
const int LP_Status_Infeasible = 1;
const int LP_Status_Unbounded = 2;
const int LP_Status_Singular = 3;
const int LP_Status_Failed = 4;

enum HSOL_INT_OPTIONS {
    INTOPT_PRINT_FLAG = 0, // 0/1 = none/do-print
    INTOPT_TRANSPOSE_FLAG, // 0/1 = none/do-transpose if possible
    INTOPT_SCALE_FLAG,     // 0/1 = none/do-scale
    INTOPT_TIGHT_FLAG,     // 0/1 = none/do-tight
    INTOPT_PERMUTE_FLAG,   // 0/1 = none/do-permute
    INTOPT_PERTURB_FLAG,   // 0/1 = none/do-perturb
    INTOPT_COUNT
};

enum HSOL_DBL_OPTIONS {
    DBLOPT_TIME_LIMIT = 0,
    DBLOPT_PRIMAL_TOL,
    DBLOPT_DUAL_TOL,
    DBLOPT_PERTURB_BASE,
    DBLOPT_PAMI_CUTOFF,
    DBLOPT_COUNT
};

enum HSOL_STR_OPTIONS {
    STROPT_PARTITION_FILE = 0, // name of row partition file
    STROPT_COUNT
};

class HModel {
public:
    /* Deakins - Classes */
    class Node{
    public:
        int data;
        Node *next;
        Node *prev;
        Node(){
            next = NULL;
            prev = NULL;
        }
    };
    class adjNode{
    public:
        int label;
        int w;
        adjNode *next;
        adjNode(){
            next = NULL;
        }
    };
    HModel();
    void setup(const char *filename);
    void setup_loadMPS(const char *filename);
    void setup_transposeLP();
    void setup_scaleMatrix();
    void setup_tightenBound();
    void setup_shuffleColumn();
    void setup_allocWorking();
    void build();

    int getNumRow() {
        return aggNumRow;
    }
    int getNumCol() {
        return aggNumCol;
    }
    int getNumTot() {
        return aggNumTot;
    }
    const HMatrix *getMatrix() {
        return &matrix;
    }
    const HFactor *getFactor() {
        return &factor;
    }
    int *getBaseIndex() {
        return &basicIndex[0];
    }
    int *getNonbasicFlag() {
        return &nonbasicFlag[0];
    }
    int *getNonbasicMove() {
        return &nonbasicMove[0];
    }
    int *getWorkIntBreak() {
        return &intBreak[0];
    }
    double *getWorkCost() {
        return &workCost[0];
    }
    double *getWorkDual() {
        return &workDual[0];
    }
    double *getWorkShift() {
        return &workShift[0];
    }
    double *getWorkLower() {
        return &workLower[0];
    }
    double *getWorkUpper() {
        return &workUpper[0];
    }
    double *getWorkRange() {
        return &workRange[0];
    }
    double *getWorkValue() {
        return &workValue[0];
    }
    double *getBaseLower() {
        return &baseLower[0];
    }
    double *getBaseUpper() {
        return &baseUpper[0];
    }
    double *getBaseValue() {
        return &baseValue[0];
    }

    void initCost(int perturb = 0);
    void initBound(int phase = 2);
    void initValue();

    void computeFactor();
    void computeDual();
    void computeDualInfeasInDual(int *dualInfeasCount);
    void computeDualInfeasInPrimal(int *dualInfeasCount);
    void correctDual(int *freeInfeasCount);
    void computePrimal();
    void computeObject(int phase = 2);

    void shiftCost(int iCol, double amount);
    void shiftBack(int iCol);
    void flipBound(int iCol);

    void updateFactor(HVector *column, HVector *row_ep, int *iRow, int *hint);
    void updateMatrix(int columnIn, int columnOut);
    void updatePivots(int columnIn, int rowOut, int sourceOut);

    void changeUpdate(int updateMethod);

    void reportPivots(int columnIn, int columnOut, double alpha);
    void reportStatus(int status);

    void printMessage(const char *message);
    void printObject();
    void printResult();
    void printProgress();
    void writePivots(const char *suffix);
    void writeMPS(const char *filename);

    // Deakins - functions
    void initStorage();
    void push(Node **headRef, int newData);
    void insertA(Node *prevNode, int newData);
    void insertB(Node **headRef, Node *nextNode, int newData);
    void append(Node **headRef, int newData);
    void appendAdj(adjNode **headRef, int newData, double newEntry);
    Node *findNode(Node **headRef, int label);
    void deleteNode(Node **headRef, Node *del);
    bool exists(Node **headRef, int data);
    int listSize(Node **headRef);
    void initEQs();
    void computeEQs();
    void lp2Graph();
    void splitColor(int s);
    void printList(Node **headRef);
    void isolate(int i);
    void equitable();
    void aggClear();   
    void aggregateA();
    void aggregateCT(int newCol);
    vector<int> getCoeff(int color);
    int getObj(int color);
    void getNewRows();
    bool discrete();


    // Solving options
    int intOption[INTOPT_COUNT];
    double dblOption[DBLOPT_COUNT];
    string strOption[STROPT_COUNT];

    // Random generator
    HRandom random;

    // The time and timer
    HTimer timer;
    double totalTime;

    // Perturbation flag
    int problemPerturbed;

    // The original models
    string modelName;

    // Solving result
    int limitUpdate;
    int countUpdate;
    int numberIteration;
    double objective;

    //private:
    // The original model
    int numCol;
    int numRow;
    int numTot;
    vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;
    double objOffset;

    /* Deakins - Objects and data structures */
    // Ints
    int masterIter;
    int vCol;
    int cCol;
    int r;
    int numParts;
    int iso = -1;   
    int prevColor; 
    int oldNumCols;
    int oldNumRows;
    int numNewRows;
    int numNewCols;
    int targ; 
    int maxColor;    
    // Storage - vectors for EPs
    vector<adjNode *> adjList;
    vector<adjNode *> aggAdjList;
    vector<double> Rhs;
    vector<Node *> C;
    vector<list<int> *> A;
    vector<double> maxcdeg;
    vector<double> mincdeg;
    vector<int> initialParts;
    vector<double> cdeg;
    vector<int> isAdj;
    vector<int> color;
    vector<bool> SCheck;
    vector<bool> isolates;
    vector<bool> linkedByIso;
    vector<int> vColor;
    vector<int> colorsToLink;
    vector<int> reps;
    vector<int> conColorReps;
    vector<int> oldColor;
    vector<int> residuals;

    // Storage - stacks
    stack<int> S;
    // Doubly linked lists
    Node *colorsAdj;
    Node *v;
    Node *c;
    Node *del;
    // Linked lists
    adjNode *w;
    // Lists
    list<int>::iterator u;
    list<int>::iterator s;
    // Storage for aggregate models
    int aggNumCol;
    int aggNumRow;
    int aggNumTot;
    vector<int> aggColIdx;
    vector<int> aggRowIdx;
    vector<int> aggAstart;
    vector<int> aggAindex;
    vector<double> aggAvalue;
    vector<double> aggColCost;
    vector<double> aggColLower;
    vector<double> aggColUpper;
    vector<double> aggRowLower;
    vector<double> aggRowUpper;
    vector<int> startingBasis;
    vector<double> startingBasicValue;
    vector<bool> prevBasicColor;
    vector<double> prevBasicValue;
    

    // Associated data of original model
    vector<int> workRowPart; // Row partition
    vector<int> intBreak;
    vector<double> dblXpert;

    // Working model
    int problemStatus;
    vector<int> basicIndex;
    vector<int> nonbasicFlag;
    vector<int> nonbasicMove;
    HMatrix matrix;
    HFactor factor;
    HVector buffer;
    HVector bufferLong;

    vector<double> workCost;
    vector<double> workDual;
    vector<double> workShift;

    vector<double> workLower;
    vector<double> workUpper;
    vector<double> workRange;
    vector<double> workValue;

    vector<double> baseLower;
    vector<double> baseUpper;
    vector<double> baseValue;

    vector<int> historyColumnIn;
    vector<int> historyColumnOut;
    vector<double> historyAlpha;
};

#endif /* HMODEL_H_ */
