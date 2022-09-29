#pragma once

#include "proxddp/core/solver-workspace.hpp"

namespace proxddp {

template <typename Scalar> struct WorkspaceFDDPTpl : WorkspaceBaseTpl<Scalar> {
  using Base = WorkspaceBaseTpl<Scalar>;
  PROXNLP_DYNAMIC_TYPEDEFS(Scalar);
  using Base::nsteps;
  using Base::q_params;
  using Base::trial_us;
  using Base::trial_xs;
  using Base::value_params;

  /// Value of `f(x_i, u_i)`
  std::vector<VectorXs> xnexts_;
  /// Feasibility gaps
  std::vector<VectorXs> feas_gaps_;
  /// State increment
  std::vector<VectorXs> dxs;
  /// Control increment
  std::vector<VectorXs> dus;
  std::vector<VectorXs> Quuks_;
  std::vector<VectorXs> ftVxx_;
  /// Buffer for KKT matrices.
  std::vector<MatrixXs> kkt_mat_bufs;
  /// Buffer for KKT system right-hand sides.
  std::vector<MatrixXs> kkt_rhs_bufs;
  /// LLT struct for each KKT system.
  std::vector<Eigen::LLT<MatrixXs>> llts_;

  explicit WorkspaceFDDPTpl(const TrajOptProblemTpl<Scalar> &problem);
};

template <typename Scalar>
WorkspaceFDDPTpl<Scalar>::WorkspaceFDDPTpl(
    const TrajOptProblemTpl<Scalar> &problem)
    : Base(problem) {

  value_params.reserve(nsteps + 1);
  q_params.reserve(nsteps);

  xnexts_.resize(nsteps + 1);
  feas_gaps_.resize(nsteps + 1);
  dxs.resize(nsteps + 1);
  dus.resize(nsteps);
  Quuks_.resize(nsteps);
  ftVxx_.resize(nsteps + 1);
  kkt_mat_bufs.resize(nsteps);
  kkt_rhs_bufs.resize(nsteps);
  llts_.reserve(nsteps);

  feas_gaps_[0] = VectorXs::Zero(problem.stages_[0]->ndx1());

  for (std::size_t i = 0; i < nsteps; i++) {
    const StageModelTpl<Scalar> &sm = *problem.stages_[i];
    const int ndx = sm.ndx1();
    const int nu = sm.nu();
    const int ndual = sm.numDual();

    value_params.emplace_back(ndx);
    q_params.emplace_back(ndx, nu, 0);

    xnexts_[i] = sm.xspace().neutral();
    feas_gaps_[i + 1] = VectorXs::Zero(sm.ndx2());
    dxs[i] = VectorXs::Zero(ndx);
    dus[i] = VectorXs::Zero(nu);

    Quuks_[i] = VectorXs::Zero(sm.nu());
    ftVxx_[i + 1] = VectorXs::Zero(sm.ndx2());
    kkt_mat_bufs[i] = MatrixXs::Zero(nu, nu);
    kkt_rhs_bufs[i] = MatrixXs::Zero(nu, ndx + 1);
    llts_.emplace_back(nu);
  }
  const StageModelTpl<Scalar> &sm = *problem.stages_.back();
  dxs[nsteps] = VectorXs::Zero(sm.ndx2());
  xnexts_[nsteps] = sm.xspace_next().neutral();
  value_params.emplace_back(sm.ndx2());

  assert(llts_.size() == nsteps);
}

} // namespace proxddp
