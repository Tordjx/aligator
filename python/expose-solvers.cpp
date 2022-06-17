#include "proxddp/python/fwd.hpp"

#include "proxddp/solver-proxddp.hpp"


namespace proxddp
{
  namespace python
  {

    void exposeSolvers()
    {
      using context::Scalar;
      using context::ShootingProblem;
      using Workspace = WorkspaceTpl<Scalar>;
      using Results = ResultsTpl<Scalar>;

      bp::class_<Workspace>(
        "Workspace", "Workspace for ProxDDP.",
        bp::init<const ShootingProblem&>(bp::args("self", "problem"))
        )
        .def_readonly("value_params", &Workspace::value_params)
        .def_readonly("q_params", &Workspace::q_params)
        .def_readonly("kkt_matrix_buffer_", &Workspace::kktMatrixFull_)
        .def_readonly("gains", &Workspace::gains_)
        .def(bp::self_ns::str(bp::self))
        ;

      bp::class_<Results>(
        "Results", "Results struct for proxDDP.", bp::init<const ShootingProblem&>()
      )
        .def_readonly("xs", &Results::xs_)
        .def_readonly("us", &Results::us_)
        .def_readonly("lams", &Results::lams_)
        .def_readonly("co_state", &Results::co_state_)
      ;

      using SolverType = SolverProxDDP<Scalar>;

      bp::enum_<MultiplierUpdateMode>("MultiplierUpdateMode", "Enum for the kind of multiplier update to use.")
        .value("NEWTON", MultiplierUpdateMode::NEWTON)
        .value("PRIMAL", MultiplierUpdateMode::PRIMAL)
        .value("PDAL", MultiplierUpdateMode::PDAL)
        ;

      bp::class_<SolverType>(
        "ProxDDP",
        bp::init< Scalar
                , Scalar, Scalar
                , Scalar, Scalar>
                ((bp::arg("self"), bp::arg("tol"),
                  bp::arg("mu_init"), bp::arg("rho_init"),
                  bp::arg("prim_alpha") = 0.1, bp::arg("prim_beta") = 0.9)
                  ))
        .def_readwrite("target_tol", &SolverType::target_tolerance, "Desired tolerance.")
        .def_readonly("mu_init",  &SolverType::mu_init,  "mu_init")
        .def_readonly("rho_init", &SolverType::rho_init, "rho_init")
        .def_readwrite("mu_factor",  &SolverType::mu_update_factor_)
        .def_readwrite("rho_factor", &SolverType::rho_update_factor_)
        .def_readwrite("multiplier_update_mode", &SolverType::mul_update_mode)
        .def("run",
             &SolverType::run,
             bp::args("self", "problem", "workspace", "results", "xs_init", "us_init"),
             "Run the solver.")
        ;
    }
    
  } // namespace python
} // namespace proxddp

