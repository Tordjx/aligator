#pragma once

#include "aligator/modelling/multibody/collision-explicit.hpp"
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/collision/distance.hpp>
#include <iostream>

namespace aligator {

template <typename Scalar>
void CollisionExplicitResidualTpl<Scalar>::evaluate(const ConstVectorRef &x,
                                                 BaseData &data) const {
  Data &d = static_cast<Data &>(data);
  pinocchio::DataTpl<Scalar> &pdata = d.pin_data_;
  Scalar pitch = x[2], yaw = x[4], posx = x[5], posy = x[6];

  // Compute rotation using Eigen
  Eigen::AngleAxis<Scalar> pitchRot(pitch, Eigen::Matrix<Scalar, 3, 1>::UnitY());
  Eigen::AngleAxis<Scalar> yawRot(yaw, Eigen::Matrix<Scalar, 3, 1>::UnitZ());
  Eigen::Quaternion<Scalar> quat = yawRot * pitchRot;

  // State vector as an Eigen vector
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> q(13);
  q << posx, posy, Scalar(0.6),
      quat.x(), quat.y(), quat.z(), quat.w(),
      0, 0, 0, 0, 0, 0;

  pinocchio::forwardKinematics(pin_model_, pdata, q);
  pinocchio::updateFramePlacements(pin_model_, pdata);

  // Computes the collision distance between pair of frames
  pinocchio::updateGeometryPlacements(pin_model_, pdata, geom_model_, d.geometry_, q);
  pinocchio::computeDistance(geom_model_, d.geometry_, frame_pair_id_);
  d.value_[0] = d.geometry_.distanceResults[frame_pair_id_].min_distance;
}

template <typename Scalar>
void CollisionExplicitResidualTpl<Scalar>::computeJacobians(const ConstVectorRef &,
                                                         BaseData &data) const {
  Data &d = static_cast<Data &>(data);
  pinocchio::DataTpl<Scalar> &pdata = d.pin_data_;

  // Calculate vector from joint to collision p1 and joint to collision p2,
  // expressed in local world aligned
  d.distance_ = d.geometry_.distanceResults[frame_pair_id_].nearest_points[0] -
                pdata.oMf[frame_id1_].translation();
  d.distance2_ = d.geometry_.distanceResults[frame_pair_id_].nearest_points[1] -
                 pdata.oMf[frame_id2_].translation();

  d.jointToP1_.setIdentity();
  d.jointToP1_.translation(d.distance_);

  d.jointToP2_.setIdentity();
  d.jointToP2_.translation(d.distance2_);

  // Get frame Jacobians
  pinocchio::computeJointJacobians(pin_model_, pdata);
  pinocchio::getFrameJacobian(pin_model_, pdata, frame_id1_,
                              pinocchio::LOCAL_WORLD_ALIGNED, d.Jcol_);

  pinocchio::getFrameJacobian(pin_model_, pdata, frame_id2_,
                              pinocchio::LOCAL_WORLD_ALIGNED, d.Jcol2_);

  // compute Jacobian at p1
  d.Jcol_ = d.jointToP1_.toActionMatrixInverse() * d.Jcol_;

  // compute Jacobian at p2
  d.Jcol2_ = d.jointToP2_.toActionMatrixInverse() * d.Jcol2_;

  // compute the residual derivatives
  // We neglect the jacobian wrt the angles of upkie

  d.Jx_.setZero();
  
  Eigen::VectorXd J_q = 
        d.geometry_.distanceResults[frame_pair_id_].normal.transpose() * 
        (d.Jcol2_.template topRows<3>() - d.Jcol_.template topRows<3>());
  d.Jx_(5,0) = J_q(0);
  d.Jx_(6,0) = J_q(1);
}

template <typename Scalar>
CollisionExplicitDataTpl<Scalar>::CollisionExplicitDataTpl(
    const CollisionExplicitResidualTpl<Scalar> &model)
    : Base(7,2, 1), pin_data_(model.pin_model_),
      geometry_(pinocchio::GeometryData(model.geom_model_)),
      Jcol_(6, model.pin_model_.nv), Jcol2_(6, model.pin_model_.nv) {
  Jcol_.setZero();
  Jcol2_.setZero();
}

} // namespace aligator
