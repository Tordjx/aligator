#pragma once

#include "aligator/modelling/multibody/multibody-wrench-cone.hpp"

#include <pinocchio/algorithm/constrained-dynamics.hpp>
#include <pinocchio/algorithm/constrained-dynamics-derivatives.hpp>

namespace aligator {

template <typename Scalar>
void MultibodyWrenchConeResidualTpl<Scalar>::evaluate(const ConstVectorRef &x,
                                                      const ConstVectorRef &u,
                                                      const ConstVectorRef &,
                                                      BaseData &data) const {
  Data &d = static_cast<Data &>(data);
  pinocchio::DataTpl<Scalar> &pdata = d.pin_data_;

  const auto q = x.head(pin_model_.nq);
  const auto v = x.tail(pin_model_.nv);

  d.tau_ = actuation_matrix_ * u;
  pinocchio::constraintDynamics(pin_model_, d.pin_data_, q, v, d.tau_,
                                constraint_models_, d.constraint_datas_,
                                d.settings);

  // Unilateral contact
  d.value_ = Acone_ * d.pin_data_.lambda_c.segment(contact_id_ * 6, 6);
}

template <typename Scalar>
void MultibodyWrenchConeResidualTpl<Scalar>::computeJacobians(
    const ConstVectorRef &x, const ConstVectorRef &, const ConstVectorRef &,
    BaseData &data) const {
  Data &d = static_cast<Data &>(data);

  const auto q = x.head(pin_model_.nq);

  pinocchio::computeConstraintDynamicsDerivatives(
      pin_model_, d.pin_data_, constraint_models_, d.constraint_datas_,
      d.settings);

  d.Jx_.leftCols(pin_model_.nv) =
      Acone_ *
      d.pin_data_.dlambda_dq.block(contact_id_ * 6, 0, 6, pin_model_.nv);
  d.Jx_.rightCols(pin_model_.nv) =
      Acone_ *
      d.pin_data_.dlambda_dv.block(contact_id_ * 6, 0, 6, pin_model_.nv);
  d.Ju_ = Acone_ *
          d.pin_data_.dlambda_dtau.block(contact_id_ * 6, 0, 6, pin_model_.nv) *
          actuation_matrix_;
}

template <typename Scalar>
MultibodyWrenchConeDataTpl<Scalar>::MultibodyWrenchConeDataTpl(
    const MultibodyWrenchConeResidualTpl<Scalar> *model)
    : Base(model->ndx1, model->nu, model->ndx2, 17),
      pin_data_(model->pin_model_), tau_(model->pin_model_.nv) {
  tau_.setZero();

  pinocchio::initConstraintDynamics(model->pin_model_, pin_data_,
                                    model->constraint_models_);
  for (auto cm = std::begin(model->constraint_models_);
       cm != std::end(model->constraint_models_); ++cm) {
    constraint_datas_.push_back(
        pinocchio::RigidConstraintDataTpl<Scalar, 0>(*cm));
  }
}

} // namespace aligator
