#pragma once

#include <proxsuite-nlp/modelling/spaces/pinocchio-groups.hpp>
#include <proxsuite-nlp/fmt-eigen.hpp>
#include <pinocchio/multibody/liegroup/special-euclidean.hpp>
#include <math.h>
#include "aligator/core/traj-opt-problem.hpp"
#include "aligator/modelling/costs/quad-state-cost.hpp"
#include "aligator/modelling/state-error.hpp"
#include "aligator/modelling/costs/sum-of-costs.hpp"
#include "aligator/modelling/dynamics/ode-abstract.hpp"
#include "aligator/modelling/dynamics/integrator-euler.hpp"

using T = double;
using VEC = proxsuite::nlp::VectorSpaceTpl<T>;

using namespace aligator;
using StateError = StateErrorResidualTpl<T>;
using QuadStateCost = QuadraticStateCostTpl<T>;
using QuadControlCost = QuadraticControlCostTpl<T>;
using context::StageModel;
using context::TrajOptProblem;
ALIGATOR_DYNAMIC_TYPEDEFS(T);

template <typename T>
struct UpkieDynamicsTpl : dynamics::ODEAbstractTpl<T> {
  using Base = dynamics::ODEAbstractTpl<T>;
  using ODEData = dynamics::ContinuousDynamicsDataTpl<T>;
  double gravity_ ; 
  double length_ ;
  UpkieDynamicsTpl(const double gravity , const double length) : Base(VEC(4), 2) {
    gravity_ = gravity;
    length_ = length;
  }


  void forward(const ConstVectorRef &x, const ConstVectorRef &u, ODEData &data) const override {
    T rdot = x[0], phidot = x[1], theta = x[2], thetadot = x[3];
    T rdotdot = u[0], phidotdot = u[1];

    data.xdot_[0] =  rdotdot;
    data.xdot_[1] = phidotdot;
    data.xdot_[2] = thetadot;
    data.xdot_[3] = std::sin(theta) * gravity_ / length_ - std::cos(theta) * rdotdot / length_;
  }

  void dForward(const ConstVectorRef &x, const ConstVectorRef &u, ODEData &data) const override {
    T theta = x[2], rdotdot = u[0];

    data.Jx_.setZero();
    data.Jx_(2,3) = 1;
    data.Jx_(3,2) =std::cos(theta) * gravity_ / length_ + std::sin(theta) * rdotdot / length_;

    data.Ju_.setZero();
    data.Ju_(0, 0) = 1;
    data.Ju_(1, 1) = 1;
    data.Ju_(3, 0) = -1 * std::cos(theta) / length_;
  }
};


