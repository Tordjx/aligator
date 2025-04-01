#pragma once

#include "aligator/modelling/multibody/collision-explicit.hpp"
#include "aligator/context.hpp"

namespace aligator {

extern template struct CollisionExplicitResidualTpl<context::Scalar>;
extern template struct CollisionExplicitDataTpl<context::Scalar>;

} // namespace aligator
