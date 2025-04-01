/// @file
/// @copyright Copyright (C) 2022-2023 LAAS-CNRS, INRIA
#ifdef ALIGATOR_WITH_PINOCCHIO
#include "aligator/fwd.hpp"
#include "aligator/python/modelling/multibody-utils.hpp"

#include "aligator/modelling/multibody/frame-placement.hpp"
#include "aligator/modelling/multibody/frame-velocity.hpp"
#include "aligator/modelling/multibody/frame-translation.hpp"
#include "aligator/modelling/multibody/frame-collision.hpp"
#include "aligator/modelling/multibody/collision-explicit.hpp"
#include "aligator/python/polymorphic-convertible.hpp"
#ifdef ALIGATOR_PINOCCHIO_V3
#include "aligator/modelling/multibody/constrained-rnea.hpp"
#endif

namespace aligator {
namespace python {
using context::ConstMatrixRef;
using context::ConstVectorRef;
using context::PinData;
using context::PinModel;
using context::StageFunction;
using context::UnaryFunction;

//
// FORWARD DECLARATIONS
//

void exposeFlyHigh();
#ifdef ALIGATOR_PINOCCHIO_V3
void exposeContactForce();
#endif
void exposeCenterOfMassFunctions();
void exposeFrameFunctions();
void exposeGravityCompensation();

//
// DEFINITIONS
//

void exposeFrameFunctions() {
  using context::Manifold;
  using context::Scalar;
  using SE3 = pinocchio::SE3Tpl<Scalar>;
  using Motion = pinocchio::MotionTpl<Scalar>;

  using FramePlacement = FramePlacementResidualTpl<Scalar>;
  using FramePlacementData = FramePlacementDataTpl<Scalar>;

  using FrameVelocity = FrameVelocityResidualTpl<Scalar>;
  using FrameVelocityData = FrameVelocityDataTpl<Scalar>;

  using FrameTranslation = FrameTranslationResidualTpl<Scalar>;
  using FrameTranslationData = FrameTranslationDataTpl<Scalar>;

  using FrameCollision = FrameCollisionResidualTpl<Scalar>;
  using FrameCollisionData = FrameCollisionDataTpl<Scalar>;

  using CollisionExplicit = CollisionExplicitResidualTpl<Scalar>;
  using CollisionExplicitData = CollisionExplicitDataTpl<Scalar>;

  using pinocchio::GeometryModel;

  if (!eigenpy::check_registration<shared_ptr<PinData>>())
    bp::register_ptr_to_python<shared_ptr<PinData>>();

  PolymorphicMultiBaseVisitor<UnaryFunction, StageFunction> unary_visitor;

  bp::class_<FramePlacement, bp::bases<UnaryFunction>>(
      "FramePlacementResidual", "Frame placement residual function.",
      bp::init<int, int, const PinModel &, const SE3 &, pinocchio::FrameIndex>(
          ("self"_a, "ndx", "nu", "model", "p_ref", "id")))
      .def(FrameAPIVisitor<FramePlacement>())
      .def(unary_visitor)
      .def("getReference", &FramePlacement::getReference, "self"_a,
           bp::return_internal_reference<>(), "Get the target frame in SE3.")
      .def("setReference", &FramePlacement::setReference, ("self"_a, "p_new"),
           "Set the target frame in SE3.");

  bp::register_ptr_to_python<shared_ptr<FramePlacementData>>();

  bp::class_<FramePlacementData, bp::bases<context::StageFunctionData>>(
      "FramePlacementData", "Data struct for FramePlacementResidual.",
      bp::no_init)
      .def_readonly("rMf", &FramePlacementData::rMf_, "Frame placement error.")
      .def_readonly("rJf", &FramePlacementData::rJf_)
      .def_readonly("fJf", &FramePlacementData::fJf_)
      .def_readonly("pin_data", &FramePlacementData::pin_data_,
                    "Pinocchio data struct.");

  bp::class_<FrameVelocity, bp::bases<UnaryFunction>>(
      "FrameVelocityResidual", "Frame velocity residual function.",
      bp::init<int, int, const PinModel &, const Motion &,
               pinocchio::FrameIndex, pinocchio::ReferenceFrame>(
          ("self"_a, "ndx", "nu", "model", "v_ref", "id", "reference_frame")))
      .def(FrameAPIVisitor<FrameVelocity>())
      .def(unary_visitor)
      .def("getReference", &FrameVelocity::getReference, "self"_a,
           bp::return_internal_reference<>(), "Get the target frame velocity.")
      .def("setReference", &FrameVelocity::setReference, ("self"_a, "v_new"),
           "Set the target frame velocity.");

  bp::register_ptr_to_python<shared_ptr<FrameVelocityData>>();

  bp::class_<FrameVelocityData, bp::bases<context::StageFunctionData>>(
      "FrameVelocityData", "Data struct for FrameVelocityResidual.",
      bp::no_init)
      .def_readonly("pin_data", &FrameVelocityData::pin_data_,
                    "Pinocchio data struct.");

  bp::class_<FrameTranslation, bp::bases<UnaryFunction>>(
      "FrameTranslationResidual", "Frame placement residual function.",
      bp::init<int, int, const PinModel &, const context::Vector3s &,
               pinocchio::FrameIndex>(
          ("self"_a, "ndx", "nu", "model", "p_ref", "id")))
      .def(FrameAPIVisitor<FrameTranslation>())
      .def(unary_visitor)
      .def("getReference", &FrameTranslation::getReference, "self"_a,
           bp::return_internal_reference<>(),
           "Get the target frame translation.")
      .def("setReference", &FrameTranslation::setReference, ("self"_a, "p_new"),
           "Set the target frame translation.");

  bp::register_ptr_to_python<shared_ptr<FrameTranslationData>>();

  bp::class_<FrameTranslationData, bp::bases<context::StageFunctionData>>(
      "FrameTranslationData", "Data struct for FrameTranslationResidual.",
      bp::no_init)
      .def_readonly("fJf", &FrameTranslationData::fJf_)
      .def_readonly("pin_data", &FrameTranslationData::pin_data_,
                    "Pinocchio data struct.");

  bp::class_<FrameCollision, bp::bases<UnaryFunction>>(
      "FrameCollisionResidual", "Frame collision residual function.",
      bp::init<int, int, const PinModel &, const GeometryModel &,
               pinocchio::PairIndex>(bp::args("self", "ndx", "nu", "model",
                                              "geom_model", "frame_pair_id")))
      .def(FrameAPIVisitor<FrameCollision>())
      .def(unary_visitor);

  bp::register_ptr_to_python<shared_ptr<FrameCollisionData>>();

  bp::class_<FrameCollisionData, bp::bases<context::StageFunctionData>>(
      "FrameCollisionData", "Data struct for FrameCollisionResidual.",
      bp::no_init)
      .def_readonly("pin_data", &FrameCollisionData::pin_data_,
                    "Pinocchio data struct.")
      .def_readonly("geom_data", &FrameCollisionData::geometry_,
                    "Geometry data struct.");


  bp::class_<CollisionExplicit, bp::bases<UnaryFunction>>(
      "CollisionExplicitResidual", "Frame collision residual function.",
      bp::init<int, int, const PinModel &, const GeometryModel &,
               pinocchio::PairIndex>(bp::args("self", "ndx", "nu", "model",
                                              "geom_model", "frame_pair_id")))
      .def(FrameAPIVisitor<CollisionExplicit>())
      .def(unary_visitor);

  bp::register_ptr_to_python<shared_ptr<CollisionExplicitData>>();

  bp::class_<CollisionExplicitData, bp::bases<context::StageFunctionData>>(
      "CollisionExplicitData", "Data struct for CollisionExplicitResidual.",
      bp::no_init)
      .def_readonly("pin_data", &CollisionExplicitData::pin_data_,
                    "Pinocchio data struct.")
      .def_readonly("geom_data", &CollisionExplicitData::geometry_,
                    "Geometry data struct.");

}

#ifdef ALIGATOR_PINOCCHIO_V3
auto underactuatedConstraintInvDyn_proxy(
    const PinModel &model, PinData &data, const ConstVectorRef &q,
    const ConstVectorRef &v, const ConstMatrixRef &actMatrix,
    const StdVectorEigenAligned<context::RCM> &constraint_models,
    StdVectorEigenAligned<context::RCD> &constraint_datas) {
  long nu = actMatrix.cols();
  int d = 0;
  for (size_t k = 0; k < constraint_models.size(); ++k) {
    d += (int)constraint_models[k].size();
  }
  context::VectorXs out(nu + d);
  underactuatedConstrainedInverseDynamics(
      model, data, q, v, actMatrix, constraint_models, constraint_datas, out);

  return bp::make_tuple((context::VectorXs)out.head(nu),
                        (context::VectorXs)out.tail(d));
}
#endif

void exposePinocchioFunctions() {
  exposeFrameFunctions();
  exposeFlyHigh();
#ifdef ALIGATOR_PINOCCHIO_V3
  exposeContactForce();
#endif
  exposeCenterOfMassFunctions();
  exposeGravityCompensation();

#ifdef ALIGATOR_PINOCCHIO_V3
  bp::def("underactuatedConstrainedInverseDynamics",
          underactuatedConstraintInvDyn_proxy,
          ("model"_a, "data", "q", "v", "actMatrix", "constraint_model",
           "constraint_data"),
          "Compute the gravity-compensating torque for a pinocchio Model under "
          "a rigid constraint.");
#endif
}
} // namespace python
} // namespace aligator

#endif