/// @file
/// @copyright Copyright (C) 2023 LAAS-CNRS, INRIA

#include "proxddp/modelling/composite-costs.hpp"

namespace aligator {

template struct CompositeCostDataTpl<context::Scalar>;

template struct QuadraticResidualCostTpl<context::Scalar>;

template struct LogResidualCostTpl<context::Scalar>;

} // namespace aligator
