/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsSeparation.h"

#include <algorithm>
#include <cassert>
#include <queue>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsModkSeparator.h"
#include "mip/HighsPathSeparator.h"
#include "mip/HighsTableauSeparator.h"
#include "mip/HighsTransformedLp.h"

HighsSeparation::HighsSeparation(const HighsMipSolver& mipsolver) {
  implBoundClock = mipsolver.timer_.clock_def("Implbound sepa", "Ibd");
  cliqueClock = mipsolver.timer_.clock_def("Clique sepa", "Clq");
  separators.emplace_back(new HighsTableauSeparator(mipsolver));
  separators.emplace_back(new HighsPathSeparator(mipsolver));
  separators.emplace_back(new HighsModkSeparator(mipsolver));
}

#if 0
class AggregationHeuristic {
  const HighsLpRelaxation& lprelaxation;
  const HighsDomain& domain;
  HighsCutPool& cutpool;
  HighsDomain& propdomain;

  const HighsMipSolver& mip;
  const HighsSolution& lpsol;
  const HighsLp& lp;

  std::vector<RowType> rowtype;
  std::vector<uint8_t> rowusable;
  std::vector<HighsInt> numcontinuous;

  std::vector<HighsInt> colubtype;
  std::vector<double> colubscale;
  std::vector<HighsInt> collbtype;
  std::vector<double> collbscale;
  std::vector<double> bounddistance;

  // vector sum object for bound substitution and retransformation
  HighsSparseVectorSum vectorsum;

  // current indices and values and rhs for the aggregated row we are looking
  // at
  std::vector<HighsInt> baseinds;
  std::vector<double> basevals;
  HighsCDouble baserhs;

  // position of continuous variables with positive bound distance
  std::vector<HighsInt> bounddistpos;

  // temporary indices and values for the row aggregating the next row
  std::vector<HighsInt> tmpinds;
  std::vector<double> tmpvals;

  // current transformed row we are looking at
  HighsCDouble rhs;
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  std::vector<double> upper;
  std::vector<double> solvals;
  std::vector<int8_t> complementation;
  HighsInt ncont;
  HighsInt nbin;
  HighsInt nint;
  HighsInt nunbndint;
  bool freevar;

  // maximal number of aggregations we allow
  static constexpr HighsInt maxAggr = 6;
  // array for the current path and its length
  HighsInt currpath[maxAggr];
  double currpathweights[maxAggr];
  HighsInt currpathlen;

  HighsInt numcuts;

 public:
  AggregationHeuristic(const HighsLpRelaxation& lprelaxation,
                       const HighsDomain& domain, HighsCutPool& cutpool,
                       HighsDomain& propdomain)
      : lprelaxation(lprelaxation),
        domain(domain),
        cutpool(cutpool),
        propdomain(propdomain),
        mip(lprelaxation.getMipSolver()),
        lpsol(lprelaxation.getLpSolver().getSolution()),
        lp(lprelaxation.getLpSolver().getLp()) {
    numcuts = 0;
  }

  void determineRowTypes() {
    rowtype.resize(lp.numRow_);
    rowusable.resize(lp.numRow_);
    for (HighsInt i = 0; i != lp.numRow_; ++i) {
      if (lp.rowLower_[i] == lp.rowUpper_[i]) {
        rowtype[i] = RowType::kEq;
        rowusable[i] = true;
        continue;
      }

      double lowerslack = kHighsInf;
      double upperslack = kHighsInf;

      if (lp.rowLower_[i] != -kHighsInf)
        lowerslack = lpsol.row_value[i] - lp.rowLower_[i];

      if (lp.rowUpper_[i] != kHighsInf)
        upperslack = lp.rowUpper_[i] - lpsol.row_value[i];

      if (lowerslack > mip.mipdata_->feastol &&
          upperslack > mip.mipdata_->feastol)
        rowtype[i] = RowType::kUnusuable;
      else if (lowerslack < upperslack)
        rowtype[i] = RowType::kGeq;
      else
        rowtype[i] = RowType::kLeq;

      rowusable[i] = rowtype[i] != RowType::kUnusuable;
    }
  }

  void detectContinuousBounds() {
    // count number of continuous variables
    numcontinuous.assign(lp.numRow_, 0);
    for (HighsInt i : mip.mipdata_->continuous_cols) {
      assert(mip.variableType(i) == HighsVarType::kContinuous);

      HighsInt start = lp.Astart_[i];
      HighsInt end = lp.Astart_[i + 1];

      for (HighsInt j = start; j != end; ++j) ++numcontinuous[lp.Aindex_[j]];
    }

    // now detect variable bound constraints with one or multiple integral
    // variables and store the sparsest one for each constinuous column Also
    // compute the bound distance to the best simple or variable bound
    colubtype.assign(lp.numCol_, -1);
    colubscale.resize(lp.numCol_);
    collbtype.assign(lp.numCol_, -1);
    collbscale.resize(lp.numCol_);
    bounddistance.assign(lp.numCol_, 0.0);

    for (HighsInt i : mip.mipdata_->continuous_cols) {
      assert(mip.variableType(i) == HighsVarType::kContinuous);

      HighsInt start = lp.Astart_[i];
      HighsInt end = lp.Astart_[i + 1];

      HighsInt lblen = kHighsIInf;
      HighsInt ublen = kHighsIInf;

      for (HighsInt j = start; j != end; ++j) {
        HighsInt row = lp.Aindex_[j];
        if (numcontinuous[row] != 1) continue;

        HighsInt rowlen = lprelaxation.getRowLen(row);
        if (rowlen != 2) continue;
        switch (rowtype[row]) {
          case RowType::kUnusuable:
            break;
          case RowType::kLeq:
            if (lp.Avalue_[j] < 0 && rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            } else if (lp.Avalue_[j] > 0 && rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            break;
          case RowType::kGeq:
            if (lp.Avalue_[j] > 0 && rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            } else if (lp.Avalue_[j] < 0 && rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            break;
          case RowType::kEq:
            if (rowlen < ublen) {
              ublen = rowlen;
              colubtype[i] = row;
              colubscale[i] = lp.Avalue_[j];
            }
            if (rowlen < lblen) {
              lblen = rowlen;
              collbtype[i] = row;
              collbscale[i] = lp.Avalue_[j];
            }
        }
      }

      double lbdist;
      double ubdist;

      assert(collbtype[i] == -1 || rowusable[collbtype[i]]);
      assert(colubtype[i] == -1 || rowusable[colubtype[i]]);
      if (collbtype[i] != -1) {
        lbdist = 0.0;
        rowusable[collbtype[i]] = false;
      } else
        lbdist = domain.colLower_[i] != -kHighsInf
                     ? lpsol.col_value[i] - domain.colLower_[i]
                     : kHighsInf;

      if (colubtype[i] != -1) {
        ubdist = 0.0;
        rowusable[colubtype[i]] = false;
      } else
        ubdist = domain.colUpper_[i] != kHighsInf
                     ? domain.colUpper_[i] - lpsol.col_value[i]
                     : kHighsInf;

      bounddistance[i] = std::min(lbdist, ubdist);
      if (bounddistance[i] < mip.mipdata_->feastol) bounddistance[i] = 0;
    }

    // todo :mark all continuous variables that have fractional integer
    // variables in their variable bound constraints then mark all rows that
    // contain fractional integer variables or marked continuous variables

    // now only count the number of continuous variables
    // not at their bound
    numcontinuous.assign(lp.numRow_, 0);
    for (HighsInt i : mip.mipdata_->continuous_cols) {
      if (bounddistance[i] == 0.0) continue;

      HighsInt start = lp.Astart_[i];
      HighsInt end = lp.Astart_[i + 1];

      for (HighsInt j = start; j != end; ++j) ++numcontinuous[lp.Aindex_[j]];
    }

    // mark each of the variables that has an equation row where it is the
    // only continuous variable by flipping its bound distance such variables
    // will be substituted in the aggregation heuristic first
  }

  void setupStartRow(HighsInt row) {
    baseinds.clear();
    basevals.clear();

    currpath[0] = row;
    currpathlen = 1;

    const HighsInt* Rindex;
    const double* Rvalue;
    HighsInt Rlen;
    lprelaxation.getRow(row, Rlen, Rindex, Rvalue);

    baseinds.insert(baseinds.end(), Rindex, Rindex + Rlen);

    if (rowtype[row] == RowType::kLeq) {
      baserhs = lp.rowUpper_[row];
      basevals.insert(basevals.end(), Rvalue, Rvalue + Rlen);
      currpathweights[0] = 1.0;
    } else {
      baserhs = -lp.rowLower_[row];
      basevals.resize(Rlen);
      currpathweights[0] = -1.0;
      for (HighsInt i = 0; i != Rlen; ++i) basevals[i] = -Rvalue[i];
    }
  }

  void computeTransformedRow(bool negate) {
    vectorsum.clear();
    bounddistpos.clear();
    inds.clear();
    vals.clear();
    upper.clear();
    solvals.clear();
    complementation.clear();

    freevar = false;

    vectorsum.setDimension(lp.numCol_);

    double basescale;

    if (negate) {
      for (HighsInt i = 0; i < currpathlen; ++i) {
        if (rowtype[currpath[i]] == RowType::kEq) continue;
        // add a slack variable for the rows in the path if they are not
        // equations
        inds.push_back(lp.numCol_ + currpath[i]);
        vals.push_back(-std::abs(currpathweights[i]));
        upper.push_back(kHighsInf);
        complementation.push_back(0);
        solvals.push_back(0);
      }
      basescale = -1.0;
    } else
      basescale = 1.0;

    rhs = basescale * baserhs;

    for (size_t j = 0; j != baseinds.size(); ++j) {
      HighsInt col = baseinds[j];
      double baseval = basescale * basevals[j];

      // add things to the vectorsum
      // store solution values and upper bounds
      // the coefficient values area added up with the sparse vector sum
      // as variable bound usage can alter things

      // we first need to complement the continuous variables
      if (mip.variableType(col) != HighsVarType::kContinuous) {
        vectorsum.add(col, baseval);
        continue;
      }

      // complement positive if possible without slack using a simple or
      // variable bound otherwise use the tightest bound and prefer simple
      // bounds adjust the right hand side accordingly
      if (bounddistance[col] > mip.mipdata_->feastol) bounddistpos.push_back(j);

      if (bounddistance[col] == kHighsInf) freevar = true;

      if (freevar) continue;

      double simplelbdist = domain.colLower_[col] != -kHighsInf
                                ? lpsol.col_value[col] - domain.colLower_[col]
                                : kHighsInf;
      double simpleubdist = domain.colUpper_[col] != kHighsInf
                                ? domain.colUpper_[col] - lpsol.col_value[col]
                                : kHighsInf;

      if (baseval < 0) {
        if (colubtype[col] != -1) {
          // use variable upper bound
          HighsInt UBrow = colubtype[col];
          double UBscale = colubscale[col];
          HighsInt UBlen;
          const HighsInt* UBindex;
          const double* UBvalue;
          double UBconst;

          if (UBscale < 0)
            UBconst = lp.rowLower_[UBrow];
          else
            UBconst = lp.rowUpper_[UBrow];

          lprelaxation.getRow(UBrow, UBlen, UBindex, UBvalue);

          HighsCDouble scale = HighsCDouble(-baseval) / UBscale;

          // the upper bound constraint adds a conmtinuous slack variable
          // with positive coefficient that we relax out

          rhs += scale * UBconst;

          for (HighsInt k = 0; k != UBlen; ++k) {
            if (UBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * UBvalue[k])) < 1e-10);
              continue;
            }

            if (mip.mipdata_->domain.colLower_[UBindex[k]] ==
                mip.mipdata_->domain.colUpper_[UBindex[k]]) {
              rhs -= mip.mipdata_->domain.colLower_[UBindex[k]] *
                     (scale * UBvalue[k]);
            } else {
              assert(mip.variableType(UBindex[k]) != HighsVarType::kContinuous);
              vectorsum.add(UBindex[k], scale * UBvalue[k]);
            }
          }
        } else if (simpleubdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use  the simple upper bound for complementation and then relax
          // the positive continuous variable as it has a positive
          // coefficient
          complementWithUpper(baseval, domain.colUpper_[col], rhs);
        } else if (simplelbdist - mip.mipdata_->feastol <= bounddistance[col]) {
          inds.push_back(col);
          complementation.push_back(1);
          if (domain.colUpper_[col] == kHighsInf)
            upper.push_back(kHighsInf);
          else
            upper.push_back(domain.colUpper_[col] - domain.colLower_[col]);
          vals.push_back(
              complementWithLower(baseval, domain.colLower_[col], rhs));
          solvals.push_back(lpsol.col_value[col] - domain.colLower_[col]);
        } else {
          assert(collbtype[col] != -1);

          // use variable lower bound
          HighsInt LBrow = collbtype[col];
          double LBscale = collbscale[col];
          HighsInt LBlen;
          const HighsInt* LBindex;
          const double* LBvalue;
          double LBconst;

          if (LBscale < 0)
            LBconst = lp.rowUpper_[LBrow];
          else
            LBconst = lp.rowLower_[LBrow];

          lprelaxation.getRow(LBrow, LBlen, LBindex, LBvalue);

          HighsCDouble scale = HighsCDouble(-baseval) / LBscale;

          rhs += scale * LBconst;

          // the lb constraint does add a slack variable with negative sign
          // as the inequality orientations do not match
          if (LBscale > 0)
            vals.push_back(-double(scale));
          else
            vals.push_back(double(scale));
          assert(vals.back() < 0);

          inds.push_back(lp.numCol_ + LBrow);
          complementation.push_back(0);
          upper.push_back(kHighsInf);
          solvals.push_back(0);

          for (HighsInt k = 0; k != LBlen; ++k) {
            if (LBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * LBvalue[k])) < 1e-10);
              continue;
            }

            if (mip.mipdata_->domain.colLower_[LBindex[k]] ==
                mip.mipdata_->domain.colUpper_[LBindex[k]]) {
              rhs -= mip.mipdata_->domain.colLower_[LBindex[k]] *
                     (scale * LBvalue[k]);
            } else {
              assert(mip.variableType(LBindex[k]) != HighsVarType::kContinuous);
              vectorsum.add(LBindex[k], scale * LBvalue[k]);
            }
          }
        }
      } else {
        assert(baseval > 0);
        if (collbtype[col] != -1) {
          HighsInt LBrow = collbtype[col];
          double LBscale = collbscale[col];
          HighsInt LBlen;
          const HighsInt* LBindex;
          const double* LBvalue;
          double LBconst;

          if (LBscale < 0)
            LBconst = lp.rowUpper_[LBrow];
          else
            LBconst = lp.rowLower_[LBrow];

          lprelaxation.getRow(LBrow, LBlen, LBindex, LBvalue);

          HighsCDouble scale = HighsCDouble(-baseval) / LBscale;

          rhs += scale * LBconst;

          assert(int(rowtype[LBrow]) * double(scale) >= 0);

          // the lower bound constraint does not add a slack variable we
          // need to remember as the inequalities match in orientation and
          // the positive slack variable is relaxed out

          for (HighsInt k = 0; k != LBlen; ++k) {
            if (LBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * LBvalue[k])) < 1e-10);
              continue;
            }

            if (mip.mipdata_->domain.colLower_[LBindex[k]] ==
                mip.mipdata_->domain.colUpper_[LBindex[k]]) {
              rhs -= mip.mipdata_->domain.colLower_[LBindex[k]] *
                     (scale * LBvalue[k]);
            } else {
              assert(mip.variableType(LBindex[k]) != HighsVarType::kContinuous);
              vectorsum.add(LBindex[k], scale * LBvalue[k]);
            }
          }
        } else if (simplelbdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use  the simple lower bound for complementation and then relax
          // the positive continuous variable as it has a positive
          // coefficient
          complementWithLower(baseval, domain.colLower_[col], rhs);
        } else if (simpleubdist - mip.mipdata_->feastol <= bounddistance[col]) {
          // use the simple upper bound and add the negative complemented
          // variable
          inds.push_back(col);
          complementation.push_back(-1);
          if (domain.colLower_[col] == -kHighsInf)
            upper.push_back(kHighsInf);
          else
            upper.push_back(domain.colUpper_[col] - domain.colLower_[col]);
          solvals.push_back(domain.colUpper_[col] - lpsol.col_value[col]);

          vals.push_back(
              complementWithUpper(baseval, domain.colUpper_[col], rhs));
        } else {
          assert(colubtype[col] != -1);
          // use variable upper bound
          HighsInt UBrow = colubtype[col];
          double UBscale = colubscale[col];
          HighsInt UBlen;
          const HighsInt* UBindex;
          const double* UBvalue;
          double UBconst;

          if (UBscale < 0)
            UBconst = lp.rowLower_[UBrow];
          else
            UBconst = lp.rowUpper_[UBrow];

          lprelaxation.getRow(UBrow, UBlen, UBindex, UBvalue);

          HighsCDouble scale = HighsCDouble(-baseval) / UBscale;
          // the ub constraint does add a slack variable as the inequality
          // orientations do not match

          if (UBscale > 0)
            vals.push_back(double(scale));
          else
            vals.push_back(-double(scale));
          assert(vals.back() < 0);

          inds.push_back(lp.numCol_ + UBrow);
          complementation.push_back(0);
          upper.push_back(kHighsInf);
          solvals.push_back(0);
          rhs += scale * UBconst;

          for (HighsInt k = 0; k != UBlen; ++k) {
            if (UBindex[k] == col)  // the column cancels out
            {
              assert(std::abs(double(baseval + scale * UBvalue[k])) < 1e-10);
              continue;
            }

            if (mip.mipdata_->domain.colLower_[UBindex[k]] ==
                mip.mipdata_->domain.colUpper_[UBindex[k]]) {
              rhs -= mip.mipdata_->domain.colLower_[UBindex[k]] *
                     (scale * UBvalue[k]);
            } else {
              assert(mip.variableType(UBindex[k]) != HighsVarType::kContinuous);
              vectorsum.add(UBindex[k], scale * UBvalue[k]);
            }
          }
        }
      }
    }

    // so far we only added continuous variables as the values of integer
    // variables could still change due to variable bound usage
    ncont = inds.size();
    nbin = 0;
    nint = 0;
    nunbndint = 0;

    // now we add the integer variables and complement them to be positive
    for (HighsInt col : vectorsum.getNonzeros()) {
      double val = vectorsum.getValue(col);

      if (std::abs(val) <= 1e-10) continue;

      if (std::abs(val) <= mip.mipdata_->feastol) {
        if (val < 0)
          rhs -= domain.colUpper_[col] * val;
        else
          rhs -= domain.colLower_[col] * val;
      }

      assert(mip.variableType(col) != HighsVarType::kContinuous);
      assert(domain.colLower_[col] != -kHighsInf ||
             domain.colUpper_[col] != kHighsInf);

      if (domain.colLower_[col] == -kHighsInf ||
          domain.colUpper_[col] == kHighsInf)
        upper.push_back(kHighsInf);
      else
        upper.push_back(
            std::floor(domain.colUpper_[col] - domain.colLower_[col] + 0.5));

      inds.push_back(col);

      if (upper.back() == 1.0)
        ++nbin;
      else if (upper.back() == kHighsInf)
        ++nunbndint;
      else
        ++nint;
      if (val < 0 && domain.colUpper_[col] != kHighsInf) {
        solvals.push_back(domain.colUpper_[col] - lpsol.col_value[col]);
        complementation.push_back(-1);
        vals.push_back(complementWithUpper(val, domain.colUpper_[col], rhs));
      } else {
        solvals.push_back(lpsol.col_value[col] - domain.colLower_[col]);
        complementation.push_back(1);
        vals.push_back(complementWithLower(val, domain.colLower_[col], rhs));
      }
    }

    // perform coefficient tightening on the transformed row
    HighsInt len = inds.size();
    HighsCDouble maxact = 0.0;
    bool unbnd = false;
    for (HighsInt i = 0; i != len; ++i) {
      if (upper[i] == kHighsInf) {
        unbnd = true;
        break;
      }
      if (vals[i] > 0) maxact += vals[i] * upper[i];
    }

    HighsCDouble maxabscoef = maxact - rhs;
    if (!unbnd && maxabscoef > mip.mipdata_->feastol) {
      HighsInt ntightened = 0;
      for (HighsInt i = 0; i != len; ++i) {
        if (inds[i] >= mip.numCol() ||
            mip.variableType(inds[i]) == HighsVarType::kContinuous)
          continue;

        if (vals[i] > maxabscoef) {
          HighsCDouble delta = vals[i] - maxabscoef;
          rhs -= delta * upper[i];
          vals[i] = double(maxabscoef);
          ++ntightened;
        }
      }
    }
  }

  void computeCut() {
    bool cutintegral;
    bool foundcut =
        generateCut(mip, upper, nbin, nint, ncont, nunbndint, solvals,
                    complementation, inds, vals, rhs, cutintegral);

    if (foundcut) {
      vectorsum.clear();
      vectorsum.setDimension(mip.numCol());

      HighsInt cutlen = inds.size();

      for (HighsInt i = 0; i != cutlen; ++i) {
        if (vals[i] == 0.0) continue;

        if (inds[i] >= mip.numCol()) {
          HighsInt row = inds[i] - mip.numCol();

          const HighsInt* Rindex;
          const double* Rvalue;
          HighsInt Rlen;
          double Rside;

          assert(rowtype[row] == RowType::kLeq || rowtype[row] == RowType::kGeq);

          lprelaxation.getRow(row, Rlen, Rindex, Rvalue);
          if (rowtype[row] == RowType::kLeq)
            Rside = lp.rowUpper_[row];
          else
            Rside = -lp.rowLower_[row];

          rhs -= vals[i] * Rside;

          double slackval = -int(rowtype[row]) * vals[i];

          for (HighsInt j = 0; j != Rlen; ++j)
            vectorsum.add(Rindex[j], slackval * Rvalue[j]);

          continue;
        }

        if (complementation[i] == 1)
          rhs += vals[i] * domain.colLower_[inds[i]];
        else if (complementation[i] == -1) {
          vals[i] = -vals[i];
          rhs += vals[i] * domain.colUpper_[inds[i]];
        }

        vectorsum.add(inds[i], vals[i]);
      }

      inds.clear();
      vals.clear();

      vectorsum.sort();

      for (HighsInt col : vectorsum.getNonzeros()) {
        double val = vectorsum.getValue(col);

        if (std::abs(val) > 1e-10) {
          vals.push_back(val);
          inds.push_back(col);
        }
      }

      double upper = double(rhs);

      cutpool.addCut(mip, inds.data(), vals.data(), inds.size(), upper,
                     cutintegral);

      ++numcuts;
    }
  }

  bool aggregateNextRow() {
    if (currpathlen == maxAggr || bounddistpos.empty()) return false;

    // sort by decreasing bound distance
    std::sort(bounddistpos.begin(), bounddistpos.end(), [&](HighsInt a, HighsInt b) {
      HighsInt col1 = baseinds[a];
      HighsInt col2 = baseinds[b];
      return bounddistance[col1] > bounddistance[col2];
    });

    HighsInt nextaggrow = -1;
    HighsCDouble nextaggscale;
    HighsInt naggbounddist = kHighsIInf;

    for (HighsInt pos : bounddistpos) {
      HighsInt col = baseinds[pos];

      HighsInt start = lp.Astart_[col];
      HighsInt end = lp.Astart_[col + 1];

      for (HighsInt j = start; j != end; ++j) {
        HighsInt row = lp.Aindex_[j];
        // if the row is marked as unusable or the inequality orientation does
        // not match we skip it
        if (!rowusable[row] ||
            basevals[pos] * lp.Avalue_[j] * int(rowtype[row]) > 0.0)
          continue;

        bool incurrentpath = false;
        for (HighsInt k = 0; k < currpathlen; ++k) {
          if (currpath[k] == row) {
            incurrentpath = true;
            break;
          }
        }

        if (incurrentpath) continue;

        bool nextaggusedbefore = currpath[0] > nextaggrow;
        bool rowusedbefore = currpath[0] > row;

        if ((nextaggusedbefore && !rowusedbefore) ||
            (nextaggusedbefore == rowusedbefore &&
             numcontinuous[row] < naggbounddist)) {
          nextaggscale = HighsCDouble(-basevals[pos]) / lp.Avalue_[j];
          nextaggrow = row;
          naggbounddist = numcontinuous[row];
          if (numcontinuous[row] == 1 && rowtype[row] == RowType::kEq) break;
        }
      }

      if (nextaggrow != -1) break;
    }

    if (nextaggrow == -1) return false;

    const HighsInt* nextRindex;
    const double* nextRvalue;
    HighsInt nextRlen;

    if (rowtype[nextaggrow] == RowType::kGeq) {
      baserhs += nextaggscale * lp.rowLower_[nextaggrow];
      assert(nextaggscale < 0);
    } else {
      assert(rowtype[nextaggrow] == RowType::kEq || nextaggscale > 0);
      baserhs += nextaggscale * lp.rowUpper_[nextaggrow];
    }

    lprelaxation.getRow(nextaggrow, nextRlen, nextRindex, nextRvalue);
    tmpinds.clear();
    tmpvals.clear();

    HighsInt a = 0;
    HighsInt b = 0;

    HighsInt baselen = baseinds.size();

    while (a != nextRlen && b != baselen) {
      if (nextRindex[a] < baseinds[b]) {
        tmpinds.push_back(nextRindex[a]);
        tmpvals.push_back(double(nextaggscale * nextRvalue[a]));
        ++a;

      } else if (baseinds[b] < nextRindex[a]) {
        tmpinds.push_back(baseinds[b]);
        tmpvals.push_back(basevals[b]);
        ++b;
      } else {
        assert(nextRindex[a] == baseinds[b]);
        double val = double(basevals[b] + nextaggscale * nextRvalue[a]);
        if (std::abs(val) > 1e-10) {
          if (std::abs(val) <= mip.mipdata_->feastol) {
            if (val < 0) {
              if (domain.colUpper_[baseinds[b]] == kHighsInf)
                return false;

              rhs -= val * domain.colUpper_[baseinds[b]];
            } else {
              if (domain.colLower_[baseinds[b]] == -kHighsInf)
                return false;

              rhs -= val * domain.colLower_[baseinds[b]];
            }
          } else {
            tmpinds.push_back(baseinds[b]);
            tmpvals.push_back(val);
          }
        }
        ++a;
        ++b;
      }
    }

    if (a != nextRlen) {
      tmpinds.insert(tmpinds.end(), nextRindex + a, nextRindex + nextRlen);
      do {
        tmpvals.push_back(double(nextaggscale * nextRvalue[a]));
        ++a;
      } while (a != nextRlen);
    }

    if (b != baselen) {
      tmpinds.insert(tmpinds.end(), baseinds.begin() + b, baseinds.end());
      tmpvals.insert(tmpvals.end(), basevals.begin() + b, basevals.end());
    }

    tmpinds.swap(baseinds);
    tmpvals.swap(basevals);
    currpath[currpathlen] = nextaggrow;
    currpathweights[currpathlen] = double(nextaggscale);
    ++currpathlen;
    return true;
  }

  void run() {
    determineRowTypes();
    detectContinuousBounds();

    for (HighsInt row = 0; row != lp.numRow_; ++row) {
      if (!rowusable[row]) continue;

      HighsInt currnumcuts = numcuts;
      setupStartRow(row);

      do {
        computeTransformedRow(false);

        if (!freevar) {
          computeCut();
          computeTransformedRow(true);
          computeCut();

          if (currnumcuts < numcuts) {
            break;
          }
        }

        if (!aggregateNextRow()) break;
      } while (true);
    }
  }
};
#endif

HighsInt HighsSeparation::separationRound(HighsDomain& propdomain,
                                          HighsLpRelaxation::Status& status) {
  const HighsSolution& sol = lp->getLpSolver().getSolution();

  HighsMipSolverData& mipdata = *lp->getMipSolver().mipdata_;

  auto propagateAndResolve = [&]() {
    if (propdomain.infeasible() || mipdata.domain.infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    propdomain.propagate();
    if (propdomain.infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    mipdata.cliquetable.cleanupFixed(mipdata.domain);
    if (mipdata.domain.infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    int numBoundChgs = (int)propdomain.getChangedCols().size();

    if (!propdomain.getChangedCols().empty()) {
      lp->setObjectiveLimit(mipdata.upper_limit);
      status = lp->resolveLp(&propdomain);

      if (!lp->scaledOptimal(status)) return -1;
    }

    return numBoundChgs;
  };

  lp->getMipSolver().timer_.start(implBoundClock);
  mipdata.implications.separateImpliedBounds(*lp, lp->getSolution().col_value,
                                             mipdata.cutpool, mipdata.feastol);
  lp->getMipSolver().timer_.stop(implBoundClock);

  HighsInt ncuts = 0;
  HighsInt numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  lp->getMipSolver().timer_.start(cliqueClock);
  mipdata.cliquetable.separateCliques(lp->getMipSolver(), sol.col_value,
                                      mipdata.cutpool, mipdata.feastol);
  lp->getMipSolver().timer_.stop(cliqueClock);

  numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  HighsTransformedLp transLp(*lp, mipdata.implications);
  if (mipdata.domain.infeasible()) {
    status = HighsLpRelaxation::Status::kInfeasible;
    return 0;
  }
  HighsLpAggregator lpAggregator(*lp);

  for (const std::unique_ptr<HighsSeparator>& separator : separators)
    separator->run(*lp, lpAggregator, transLp, mipdata.cutpool);

  numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  mipdata.cutpool.separate(sol.col_value, propdomain, cutset, mipdata.feastol);

  if (cutset.numCuts() > 0) {
    ncuts += cutset.numCuts();
    lp->addCuts(cutset);
    status = lp->resolveLp(&propdomain);
    lp->performAging();
  }

  return ncuts;
}

void HighsSeparation::separate(HighsDomain& propdomain) {
  HighsLpRelaxation::Status status = lp->getStatus();
  const HighsMipSolver& mipsolver = lp->getMipSolver();

  if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty()) {
    // double firstobj = lp->getObjective();
    double firstobj = mipsolver.mipdata_->rootlpsolobj;

    while (lp->getObjective() < mipsolver.mipdata_->upper_limit) {
      double lastobj = lp->getObjective();

      size_t nlpiters = -lp->getNumLpIterations();
      HighsInt ncuts = separationRound(propdomain, status);
      nlpiters += lp->getNumLpIterations();
      mipsolver.mipdata_->sepa_lp_iterations += nlpiters;
      mipsolver.mipdata_->total_lp_iterations += nlpiters;
      // printf("separated %" HIGHSINT_FORMAT " cuts\n", ncuts);

      // printf(
      //     "separation round %" HIGHSINT_FORMAT " at node %" HIGHSINT_FORMAT "
      //     added %" HIGHSINT_FORMAT " cuts objective changed " "from %g to %g,
      //     first obj is %g\n", nrounds, (HighsInt)nnodes, ncuts, lastobj,
      //     lp->getObjective(), firstobj);
      if (ncuts == 0 || !lp->scaledOptimal(status) ||
          lp->getFractionalIntegers().empty())
        break;

      // if the objective improved considerably we continue
      if ((lp->getObjective() - firstobj) <=
          std::max((lastobj - firstobj), mipsolver.mipdata_->feastol) * 1.01)
        break;
    }

    // printf("done separating\n");
  } else {
    // printf("no separation, just aging. status: %" HIGHSINT_FORMAT "\n",
    // (HighsInt)status);
    lp->performAging(false);
    mipsolver.mipdata_->cutpool.performAging();
  }
}
