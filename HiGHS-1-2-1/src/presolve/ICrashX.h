/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_ICRASHX_H_
#define PRESOLVE_ICRASHX_H_

#include <iostream>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"

bool callCrossover(const HighsLp& lp, const HighsOptions& options,
                   HighsSolution& solution, HighsBasis& highs_basis, HighsInfo& highs_info);
bool callPrimalCrossover(const HighsLp& lp, const HighsOptions& options,
                   HighsSolution& solution, HighsBasis& highs_basis, HighsInfo& highs_info);            

#endif
