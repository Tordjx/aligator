/// @file
/// @copyright Copyright (C) 2022 LAAS-CNRS, INRIA
#pragma once

#include "proxddp/core/function-abstract.hpp"

#include <proxnlp/modelling/spaces/vector-space.hpp>

#include "proxddp/core/cost-abstract.hpp"
#include "proxddp/core/dynamics.hpp"
#include "proxddp/core/constraint.hpp"
#include "proxddp/utils/exceptions.hpp"

#include "proxddp/core/clone.hpp"

namespace proxddp {

/** @brief    A stage in the control problem.
 *
 *  @details  Each stage containts cost functions, dynamical
 *            and constraint models. These objects are hold
 *            through smart pointers to leverage dynamic
 *            polymorphism.
 */
template <typename _Scalar>
struct StageModelTpl : Cloneable<StageModelTpl<_Scalar>> {
public:
  using Scalar = _Scalar;
  PROXNLP_DYNAMIC_TYPEDEFS(Scalar);

  using Manifold = ManifoldAbstractTpl<Scalar>;
  using ManifoldPtr = shared_ptr<Manifold>;
  using Dynamics = DynamicsModelTpl<Scalar>;
  using DynamicsPtr = shared_ptr<Dynamics>;
  using FunctionPtr = shared_ptr<StageFunctionTpl<Scalar>>;
  using ConstraintSetPtr = shared_ptr<ConstraintSetBase<Scalar>>;
  using Constraint = StageConstraintTpl<Scalar>;
  using Cost = CostAbstractTpl<Scalar>;
  using Data = StageDataTpl<Scalar>;
  using VectorSpace = proxnlp::VectorSpaceTpl<Scalar>;

  /// State space for the current state \f$x_k\f$.
  ManifoldPtr xspace_;
  /// State space for the next state \f$x_{k+1}\f$.
  ManifoldPtr xspace_next_;
  /// Control vector space -- by default, a simple Euclidean space.
  ManifoldPtr uspace_;
  /// Stage cost function.
  shared_ptr<Cost> cost_;
  /// Constraint manager.
  ConstraintStackTpl<Scalar> constraints_;

  const Manifold &xspace() const { return *xspace_; }
  const Manifold &uspace() const { return *uspace_; }
  const Manifold &xspace_next() const { return *xspace_next_; }
  virtual const Dynamics &dyn_model() const {
    assert(constraints_.numConstraints() > 0);
    auto dyn_ptr = std::static_pointer_cast<Dynamics>(constraints_[0].func);
    if (dyn_ptr == 0) {
      PROXDDP_RUNTIME_ERROR(
          "First element of constraints_ array should be a DynamicsModel.");
    } else {
      return *dyn_ptr;
    }
  }

  const Constraint &getConstraint(std::size_t j) const {
    if (j >= constraints_.numConstraints())
      PROXDDP_RUNTIME_ERROR("Maximum index exceeded.");
    return constraints_[j];
  }

  virtual const Cost &cost() const { return *cost_; }

  int nx1() const { return xspace_->nx(); }
  int ndx1() const { return xspace_->ndx(); }
  int nu() const { return uspace_->ndx(); }
  int nx2() const { return xspace_next_->nx(); }
  int ndx2() const { return xspace_next_->ndx(); }

  /// Number of constraints (constraint objects).
  std::size_t numConstraints() const { return constraints_.numConstraints(); }

  /// Number of primal optimization variables.
  int numPrimal() const;
  /// Number of dual variables, i.e. Lagrange multipliers.
  int numDual() const;

  /// Constructor assumes the control space is a Euclidean space of
  /// dimension \p nu.
  StageModelTpl(shared_ptr<Cost> cost, DynamicsPtr dyn_model);
  virtual ~StageModelTpl() = default;

  /// @brief    Add a constraint to the stage.
  template <typename T> void addConstraint(T &&cstr);

  /// @copybrief  addConstraint().
  /// @details    Adds a constraint by allocating a new StageConstraintTpl.
  void addConstraint(FunctionPtr func, ConstraintSetPtr cstr_set);

  /* Evaluate costs, constraints, ... */

  /// @brief    Evaluate all the functions (cost, dynamics, constraints) at this
  /// node.
  virtual void evaluate(const ConstVectorRef &x, const ConstVectorRef &u,
                        const ConstVectorRef &y, Data &data) const;

  /// @brief    Compute the derivatives of the StageModelTpl.
  virtual void computeDerivatives(const ConstVectorRef &x,
                                  const ConstVectorRef &u,
                                  const ConstVectorRef &y, Data &data) const;

  /// @brief    Create a Data object.
  virtual shared_ptr<Data> createData() const;

  template <typename S>
  friend std::ostream &operator<<(std::ostream &oss,
                                  const StageModelTpl<S> &stage);

protected:
  StageModelTpl(ManifoldPtr space, const int nu);
};

/// @brief    Data struct for stage models StageModelTpl.
template <typename _Scalar>
struct StageDataTpl : public Cloneable<StageDataTpl<_Scalar>> {
  using Scalar = _Scalar;
  PROXNLP_DYNAMIC_TYPEDEFS(Scalar);

  using StageModel = StageModelTpl<Scalar>;
  using CostDataAbstract = CostDataAbstractTpl<Scalar>;
  using FunctionData = FunctionDataTpl<Scalar>;
  using DynamicsData = DynamicsDataTpl<Scalar>;

  /// Data structs for the functions involved in the constraints.
  std::vector<shared_ptr<FunctionData>> constraint_data;
  /// Data for the running costs.
  shared_ptr<CostDataAbstract> cost_data;

  /// @brief    Constructor.
  ///
  /// @details  The constructor initializes or fills in the data members using
  /// move semantics.
  explicit StageDataTpl(const StageModel &stage_model);

  virtual ~StageDataTpl() = default;

  DynamicsData &dyn_data() {
    return static_cast<DynamicsData &>(*constraint_data[0]);
  }

  const DynamicsData &dyn_data() const {
    return static_cast<const DynamicsData &>(*constraint_data[0]);
  }

  /// @brief Check data integrity.
  virtual void checkData() {
    const char msg[] = "StageData integrity check failed.";
    if (constraint_data.size() == 0) {
      PROXDDP_RUNTIME_ERROR(fmt::format("{} (constraint_data empty)", msg));
    }
    if (cost_data == 0) {
      PROXDDP_RUNTIME_ERROR(fmt::format("{} (cost_data is nullptr)", msg));
    }
    shared_ptr<const DynamicsData> dd =
        std::dynamic_pointer_cast<const DynamicsData>(constraint_data[0]);
    if (dd == nullptr) {
      PROXDDP_RUNTIME_ERROR(
          fmt::format("{} (constraint_data[0] should be dynamics data)", msg));
    }
  }

protected:
  StageDataTpl() = default;
};

} // namespace proxddp

#include "proxddp/core/stage-model.hxx"
