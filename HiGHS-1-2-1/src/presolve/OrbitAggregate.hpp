#ifndef ORBIT_AGGREGATE_H_
#define ORBIT_AGGREGATE_H_

#include "orbital_partition.hpp"
#include "lp_data/HighsLp.h"
#include <set>
#include <map>
#include <limits>

class OrbitAggregate{
public:
    void setup(orbital_partition& orbits, HighsLp& lp, 
                   std::vector<HighsInt>& parent_links, std::vector<HighsInt>& child_links);
    void updateMasterLpAndEp(orbital_partition& op, int _nC, int _nR,
                            int _nnz, std::vector<int>& As, std::vector<int>& Ai,
                            std::vector<double>& Av, std::vector<double>& rL, std::vector<double>& rU);
    void clear();
    void findDimensions();
    void buildResidualRows();
    void pairCutWithIndex();
    void pairOrderConstraintWithIndex();
    void aggregateColBounds();
    void aggregateRowBounds();
    void aggregateAMatrix();
    void addCutsToAggregate();
    void addOrderConstraints();
    void aggregateCostVector();
    void aggregate();
    void buildBasis();
    int getNumCol();
    int getNumRowAfterCuts();
    int getNumRowAfterOrderConstraints();
    int getNumRow();
    HighsLp& getAggLp() { return agg_lp_; }
    std::vector<double>& getColUpper();
    std::vector<double>& getColLower();
    std::vector<double>& getColCost();
    std::vector<double>& getRowUpper();
    std::vector<double>& getRowLower();
    std::vector<double>& getAvalue();
    std::vector<int>& getAindex();
    std::vector<int>& getAstart();

    // Orbital partition container
    orbital_partition orbits_;
    // Lp containers
    HighsLp original_lp_;
    HighsLp agg_lp_;
    // Residual col/row containers
    std::vector<HighsInt> parent_links_;
    std::vector<HighsInt> child_links_;
    std::vector<HighsInt> parent_row_;
    std::vector<HighsInt> child_row_;
    std::vector<HighsInt> is_parent_;
    std::vector<HighsInt> is_child_;
    std::vector<HighsInt> parent_scale_;
    std::vector<HighsInt> residual_row_;
    // Basis container
    HighsBasis opt_basis_;
    HighsBasis basis_;

};

#endif 