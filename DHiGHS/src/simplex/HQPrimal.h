/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HQPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HQPRIMAL_H_
#define SIMPLEX_HQPRIMAL_H_

#include <utility>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"
#include "presolve/hopcroftKarp.h"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Dense"

using std::pair;


/**
 * @brief Phase 2 primal simplex solver for HiGHS
 *
 * Not an efficient primal simplex solver: just a way of tidying up
 * dual infeasibilities when dual optimality (primal feasibility) has
 * been acheived with the dual simplex method
 */
class HQPrimal {
 public:
  HQPrimal(HighsModelObject& model_object) : workHMO(model_object) {}
  /**
   * @brief Solve a model instance
   */
  HighsStatus solve();

  /**
   * @brief Perform Phase 2 primal simplex iterations
   */
  void solvePhase2();
  void solvePhase3();
  void solvePhaseSwap();

 private:
  void primalRebuild();
  void primalChooseColumn();
  void unfold();
  void swapDegenerate();
  void primalChooseR();
  void primalChooseDegenerateSlack();
  void computeReduceResiduals(int i);
  void buildDependencyGraph(int slack);
  void buildDependencyMatrix();
  void buildSlackRSubMatrix();
  void buildGeneralQR();
  void buildReduceLUFactor();
  void buildZeroMatching();
  void buildMaximumMatching();
  void matchingHeuristic();
  void buildNewBasis();
  void buildNewQRBasis();
  void buildNewLUBasis();
  void changeDegenSlack();
  void appendLeftOverResiduals();
  void buildInitialFactor();
  void updateSimplexLp();
  void updateHighsBasis();
  void updateSolver();
  void updateAMatrix();
  void countNonBoundedVars();
  void primalChooseRow();
  void primalUpdate();
  int isReadyForHighs();

  void phase1ComputeDual();
  void phase1ChooseColumn();
  void phase1ChooseColumnArtificial();
  void phase1ChooseRow();
  void phase1Update();

  void devexReset();
  void devexUpdate();

  void buildTableau();

  /**
   * @brief Pass the data for the iteration analysis, report and rebuild report
   */
  void iterationAnalysisData();

  /**
   * @brief Perform the iteration analysis
   */
  void iterationAnalysis();

  /**
   * @brief Single line report after rebuild
   */
  void reportRebuild(const int rebuild_invert_hint=-1);

  // Model pointer
  HighsModelObject& workHMO;
  hopcroftKarp hK;

  int solver_num_col;
  int originalNumCol;
  int solver_num_row;
  int originalNumRow;
  int solver_num_tot;
  int solver_num_between_bounds;
  int solver_num_at_bounds;
  int readyForHighs = 0;
  HighsSimplexAnalysis* analysis;

  bool no_free_columns;

  int isPrimalPhase1;

  int solvePhase;
  // Pivot related
  int numDrop;
  int invertHint;
  int columnIn;
  int rowOut;
  int columnOut;
  int phase1OutBnd;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  //  double alphaRow;
  double numericalTrouble;
  int num_flip_since_rebuild;

  // Primal phase 1 tools
  vector<pair<double, int> > ph1SorterR;
  vector<pair<double, int> > ph1SorterT;

  // Devex weight
  int num_devex_iterations;
  int num_bad_devex_weight;
  vector<double> devex_weight;
  vector<int> devex_index;

  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;

  // Eigen Stuff for factorization
  vector<int> depStart;
  vector<int> depInd;
  vector<int> depMatStart;
  vector<int> depMatIndex;
  vector<int> depMatValue;
  vector<Eigen::Triplet<double> > spTriplets;
  Eigen::SparseMatrix<double> spDepMat;
  Eigen::SparseMatrix<double> spSubMat;
  // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > LU;
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > QR;
  Eigen::FullPivLU<Eigen::SparseMatrix<double> > LU;
  const int* spSubMatQRPerm;
  std::vector<int> spSubMatRowPerm;
  std::vector<int> spSubMatRowUnperm;
  std::vector<int> spSubMatColPerm;
  std::vector<int> spSubMatColUnperm;
  std::vector<int> rColSwapped;
  std::vector<std::vector<double> > printBinvAr;
  int rowPerms = -1;
  int colPerms = -1;
  int spSubMatRank;
  int spSubMatCols;
  matching swap;
  // New method using highs LU factorization
  int rReduceNumCol = 0;
  int rReduceNumRow = 0;
  int rReduceNnz = 0;
  std::vector<int> rReduceAstart;
  std::vector<int> rZeroStart;
  std::vector<int> rReduceAindex;
  std::vector<int> rZeroIndex;
  std::vector<double> rReduceAvalue;
  std::vector<int> rSwapBasis;
  std::vector<int> rReduceColPerm;
  std::vector<int> degenRowPair;
  std::vector<int> rColPair;
  std::vector<int> zeros;
  int numPairs = 0;
  HFactor rReduceFactor;
  int rReduceRankDeficiency;
  int rReduceRank;
  HMatrix rReduceMatrix;

  // Artificial variable tracker
  vector<bool> pivotArtificial;

  double row_epDensity;
  double columnDensity;
  bool degeneratePivot;
  std::vector<bool> skipRow;
  std::vector<bool> rColPaired;
  std::vector<bool> slackPaired;
  std::vector<bool> rInNeighborhood;
  std::vector<int> rColPairing;
  std::vector<int> slackPairing;

  // Containers for variable swapping in degenerate submatrix
  std::vector<int> rSwapped;
  std::vector<int> sSwapped;
};

#endif /* SIMPLEX_HQPRIMAL_H_ */
