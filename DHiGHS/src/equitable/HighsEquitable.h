#ifndef HIGHS_EQUITABLE_H
#define HIGHS_EQUITABLE_H

#include "HighsLp.h"
#include "HighsSimpleDec.h"
#include "HighsTimer.h"
#include "saucy-equitable.h"
#include "util.h"
#include "amorph.h"
#include "platform.h"

#include <string>
#include <vector>
#include <algorithm>
#include <stack>
#include <set>
#include <list>
#include <map>
#include <tuple>
#include <numeric>
#include <functional>
#include <forward_list>
#include <signal.h>
using namespace std;

class HighsEquitable {
public:
	// Setup for equitable ptn
	HighsEquitable(const HighsLp& lp);
    eq_part* refine();
    void lp2Graph();
    void doSaucyEquitable();
    int init_fixadj1(int n, int* adj);
    void init_fixadj2(int n, int e, int* adj);
    static int on_automorphism( int n, const int *gamma, int k, int *support, void *arg );
    static void amorph_print_automorphism(
    int n, const int *gamma, int nsupp, const int *support,
    struct amorph_graph *g, char *marks );

    // Saucy
    struct saucy *s;
    struct eq_part *partitions;
    struct saucy_stats stats;
	struct amorph_graph *g;
    const static int quiet_mode = 0;
    const static sig_atomic_t timeout_flag = 0;
    int timeout = 0;
    static char *marks;
    
	/// Original LP info
    int nRows;
    int nCols; 
    int nTot;
    static int nTotal;
    int nnz;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper; 
    vector<double> rowLower;
    vector<double> rowUpper;
    vector<double> Avalue; 
    vector<int> Aindex; 
    vector<int> Astart;
    vector<double> AvaluePos;
    vector<double> ARvalue;
    vector<int> ARindex;
    vector<int> ARstart;
    vector<double> ARvaluePos;
    vector<int> AindexP;
    vector<string> rowNames;
    vector<string> colNames;
};

#endif