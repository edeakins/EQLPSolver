#include "OrbitAggregate.hpp"

void OrbitAggregate::setup(orbital_partition& orbits, HighsLp& original_lp,
                               std::vector<HighsInt>& parent_links, std::vector<HighsInt>& child_links){
    orbits_ = orbits;
    original_lp_ = original_lp;
    parent_links_ = parent_links;
    child_links_ = child_links;
}

void OrbitAggregate::aggregate(){
    findDimensions();
    buildResidualRows();
    aggregateColBounds();
    aggregateRowBounds();
    aggregateAMatrix();
    aggregateCostVector();
}

void OrbitAggregate::findDimensions(){
    int i_part, el, el_idx;
    agg_lp_.num_col_ = agg_lp_.num_row_ = 0;
    for (i_part = 0; i_part < orbits_.orbit_start.size() - 1; ++i_part){
        el_idx = orbits_.orbit_start.at(i_part);
        el = orbits_.element.at(el_idx);
        if (el < original_lp_.num_col_)
            agg_lp_.num_col_++;
        else
            agg_lp_.num_row_++;;
    }
    agg_lp_.num_aggregate_cols_ = agg_lp_.num_col_;
    agg_lp_.num_aggregate_rows_ = agg_lp_.num_row_;
    agg_lp_.num_col_ += parent_links_.size();
    agg_lp_.num_row_ += parent_links_.size();
    agg_lp_.a_matrix_.num_col_ = agg_lp_.num_col_;
    agg_lp_.a_matrix_.num_row_ = agg_lp_.num_row_;
    agg_lp_.num_residual_cols_ = parent_links_.size();
    agg_lp_.num_residual_rows_ = parent_links_.size();
    agg_lp_.col_lower_.resize(agg_lp_.num_col_);
    agg_lp_.col_upper_.resize(agg_lp_.num_col_);
    agg_lp_.col_cost_.resize(agg_lp_.num_col_); 
    agg_lp_.row_lower_.resize(agg_lp_.num_row_);
    agg_lp_.row_upper_.resize(agg_lp_.num_row_);
}

void OrbitAggregate::buildResidualRows(){
    if (!agg_lp_.num_residual_cols_) return;
    HighsInt res_row, par, chi;
    parent_row_.resize(agg_lp_.num_aggregate_cols_ + 1);
    is_parent_.resize(agg_lp_.num_aggregate_cols_);
    child_row_.resize(agg_lp_.num_aggregate_cols_);
    is_child_.resize(agg_lp_.num_aggregate_cols_);
    parent_scale_.clear();
    std::vector<HighsInt> parent_freq(agg_lp_.num_aggregate_cols_);
    HighsInt min_parent = kHighsIInf;
    for (HighsInt i = 0; i < parent_links_.size(); ++i){
        HighsInt i_col = parent_links_.at(i);
        if (!parent_freq.at(i_col)++) is_parent_.at(i_col) = 1;
        if (i_col < min_parent) min_parent = i_col;
    }
    parent_row_.at(min_parent) = agg_lp_.num_aggregate_rows_;
    for (HighsInt i_col = min_parent; i_col < agg_lp_.num_aggregate_cols_; ++i_col)
        parent_row_.at(i_col + 1) = parent_row_.at(i_col) + parent_freq.at(i_col);
    for (res_row = 0; res_row < agg_lp_.num_residual_rows_; ++res_row){
        par = parent_links_.at(res_row);
        chi = child_links_.at(res_row);
        is_child_.at(chi) = 1;
        child_row_.at(chi) = agg_lp_.num_aggregate_rows_ + res_row;
        double par_orb_size = orbits_.orbit_start.at(par + 1) -
                          orbits_.orbit_start.at(par);
        double chi_orb_size = orbits_.orbit_start.at(chi + 1) - 
                          orbits_.orbit_start.at(chi);
        parent_scale_.push_back((double)chi_orb_size/par_orb_size);
        residual_row_.push_back(agg_lp_.num_aggregate_rows_ + res_row);
    }
}

void OrbitAggregate::aggregateColBounds(){
    HighsInt i_col;
    for (i_col = 0; i_col < agg_lp_.num_aggregate_cols_; ++i_col){
        HighsInt el_idx = orbits_.orbit_start.at(i_col);
        HighsInt el = orbits_.element.at(el_idx);
        HighsInt orbit_size = orbits_.orbit_start.at(i_col + 1) - el_idx;
        agg_lp_.col_lower_.at(i_col) = original_lp_.col_lower_.at(el) * orbit_size;
        agg_lp_.col_upper_.at(i_col) = original_lp_.col_upper_.at(el) * orbit_size;
    }
}

void OrbitAggregate::aggregateRowBounds(){
    HighsInt i_row;
    for (i_row = 0; i_row < agg_lp_.num_aggregate_rows_; ++i_row){
        HighsInt el_idx = orbits_.orbit_start.at(i_row + agg_lp_.num_aggregate_cols_);
        HighsInt el = orbits_.element.at(el_idx) - original_lp_.num_col_;
        HighsInt orbit_size = (orbits_.orbit_start.at(i_row + 
                                                      agg_lp_.num_aggregate_cols_
                                                       + 1) - el_idx);
        if (original_lp_.row_lower_.at(el) == kHighsInf)
            agg_lp_.row_lower_.at(i_row) = kHighsInf;
        else 
            agg_lp_.row_lower_.at(i_row) = (original_lp_.row_lower_.at(el)
                                            * orbit_size);
        if (original_lp_.row_upper_.at(el) == kHighsInf)
            agg_lp_.row_upper_.at(i_row) = kHighsInf;
        else 
            agg_lp_.row_upper_.at(i_row) = (original_lp_.row_upper_.at(el)
                                            * orbit_size);
    }
}

void OrbitAggregate::aggregateAMatrix(){
    HighsInt i_col, i_mat, r_col;
    for (i_col = 0; i_col < agg_lp_.num_aggregate_cols_; ++i_col){
        HighsInt el_idx = orbits_.orbit_start.at(i_col);
        HighsInt el = orbits_.element.at(el_idx);
        std::vector<double> coeff(agg_lp_.num_aggregate_rows_, 0);
        for (i_mat = original_lp_.a_matrix_.start_.at(el); 
             i_mat < original_lp_.a_matrix_.start_.at(el + 1); ++i_mat){
            HighsInt a_idx = original_lp_.a_matrix_.index_.at(i_mat);
            HighsInt row_idx = orbits_.orbit.at(a_idx + original_lp_.num_col_) - agg_lp_.num_aggregate_cols_;
            coeff.at(row_idx) += original_lp_.a_matrix_.value_.at(i_mat);
        }
        for (i_mat = 0; i_mat < coeff.size(); ++i_mat){
            if (coeff.at(i_mat)){
                agg_lp_.a_matrix_.value_.push_back(coeff.at(i_mat));
                agg_lp_.a_matrix_.index_.push_back(i_mat);
            }
        }
        if (agg_lp_.num_residual_cols_ && is_parent_.at(i_col)){
            for (HighsInt agg_row = parent_row_.at(i_col); agg_row < parent_row_.at(i_col + 1); ++agg_row){
                HighsInt offset = agg_row - agg_lp_.num_aggregate_rows_;
                double scale = parent_scale_.at(offset); 
                agg_lp_.a_matrix_.value_.push_back(scale);
                agg_lp_.a_matrix_.index_.push_back(agg_row);
            }
        }
        if (agg_lp_.num_residual_cols_ && is_child_.at(i_col)){
            double coeff = orbits_.orbit_start.at(i_col + 1) -
                orbits_.orbit_start.at(i_col);
            agg_lp_.a_matrix_.value_.push_back(-1.0);
            agg_lp_.a_matrix_.index_.push_back(child_row_.at(i_col));
        }
        agg_lp_.a_matrix_.start_.push_back(agg_lp_.a_matrix_.value_.size());
    }
    for (r_col = 0; r_col < agg_lp_.num_residual_cols_; ++r_col){
        agg_lp_.a_matrix_.value_.push_back(-1.0);
        agg_lp_.a_matrix_.index_.push_back(residual_row_.at(r_col));
        agg_lp_.a_matrix_.start_.push_back(agg_lp_.a_matrix_.value_.size());
    }
}

void OrbitAggregate::aggregateCostVector(){
    HighsInt i_col;
    for (i_col = 0; i_col < agg_lp_.num_aggregate_cols_; ++i_col){
        HighsInt el_idx = orbits_.orbit_start.at(i_col);
        HighsInt el = orbits_.element.at(el_idx);
        agg_lp_.col_cost_.at(i_col) = original_lp_.col_cost_.at(el);
    }
}