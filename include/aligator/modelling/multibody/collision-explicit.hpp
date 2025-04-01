#pragma once

#include "aligator/core/unary-function.hpp"
#include "./fwd.hpp"

#include <pinocchio/multibody/geometry.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/frame.hpp>
#include <pinocchio/algorithm/geometry.hpp>

#include <proxsuite-nlp/third-party/polymorphic_cxx14.hpp>

namespace aligator {

template <typename Scalar> struct CollisionExplicitDataTpl;

template <typename _Scalar>
struct CollisionExplicitResidualTpl : UnaryFunctionTpl<_Scalar>, frame_api {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using Scalar = _Scalar;
  ALIGATOR_DYNAMIC_TYPEDEFS(Scalar);
  ALIGATOR_UNARY_FUNCTION_INTERFACE(Scalar);
  using BaseData = typename Base::Data;
  using Model = pinocchio::ModelTpl<Scalar>;
  using ManifoldPtr = xyz::polymorphic<ManifoldAbstractTpl<Scalar>>;
  using SE3 = pinocchio::SE3Tpl<Scalar>;
  using Data = CollisionExplicitDataTpl<Scalar>;
  using GeometryModel = pinocchio::GeometryModel;

  Model pin_model_;
  GeometryModel geom_model_;

  CollisionExplicitResidualTpl(const int ndx, const int nu, const Model &model,
                            const GeometryModel &geom_model,
                            const pinocchio::PairIndex frame_pair_id)
      : Base(7, 2, 1), pin_model_(model), geom_model_(geom_model),
        frame_pair_id_(frame_pair_id) {
    frame_id1_ =
        geom_model
            .geometryObjects[geom_model.collisionPairs[frame_pair_id_].first]
            .parentFrame;
    frame_id2_ =
        geom_model
            .geometryObjects[geom_model.collisionPairs[frame_pair_id_].second]
            .parentFrame;
  }

  void evaluate(const ConstVectorRef &x, BaseData &data) const;

  void computeJacobians(const ConstVectorRef &x, BaseData &data) const;

  shared_ptr<BaseData> createData() const {
    return allocate_shared_eigen_aligned<Data>(*this);
  }

protected:
  pinocchio::PairIndex frame_pair_id_;
  pinocchio::FrameIndex frame_id1_;
  pinocchio::FrameIndex frame_id2_;
};

template <typename Scalar>
struct CollisionExplicitDataTpl : StageFunctionDataTpl<Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using Base = StageFunctionDataTpl<Scalar>;
  using PinData = pinocchio::DataTpl<Scalar>;
  using PinGeom = pinocchio::GeometryData;
  using pinPlacement = pinocchio::SE3;

  /// Pinocchio data object
  PinData pin_data_;
  /// Pinocchio geometry object
  pinocchio::GeometryData geometry_;
  /// Jacobian of the collision point
  typename math_types<Scalar>::Matrix6Xs Jcol_;
  typename math_types<Scalar>::Matrix6Xs Jcol2_;
  /// Placement of collision point to joint
  pinPlacement jointToP1_;
  pinPlacement jointToP2_;
  /// Distance from nearest point to joint for each collision frame
  typename math_types<Scalar>::Vector3s distance_;
  typename math_types<Scalar>::Vector3s distance2_;

  CollisionExplicitDataTpl(const CollisionExplicitResidualTpl<Scalar> &model);
};

} // namespace aligator

#ifdef ALIGATOR_ENABLE_TEMPLATE_INSTANTIATION
#include "aligator/modelling/multibody/collision-explicit.txx"
#endif
