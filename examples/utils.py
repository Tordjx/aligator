"""
Common utilities for examples.
"""
import numpy as np
import pinocchio as pin
import hppfcl as fcl
import example_robot_data as erd
import tap

import matplotlib.pyplot as plt
from typing import Literal, List


plt.rcParams["lines.linewidth"] = 1.0
plt.rcParams["lines.markersize"] = 5

integrator_choices = Literal["euler", "semieuler", "midpoint", "rk2"]


class ArgsBase(tap.Tap):
    display: bool = False  # Display the trajectory using meshcat
    record: bool = False  # record video
    plot: bool = False
    integrator: integrator_choices = "semieuler"
    """Numerical integrator to use"""


def plot_se2_pose(q: np.ndarray, ax: plt.Axes, alpha=0.5, fc="tab:blue"):
    from matplotlib import transforms

    w = 1.0
    h = 0.4
    center = (q[0] - 0.5 * w, q[1] - 0.5 * h)
    rect = plt.Rectangle(center, w, h, fc=fc, alpha=alpha)
    theta = np.arctan2(q[3], q[2])
    transform_ = transforms.Affine2D().rotate_around(*q[:2], -theta) + ax.transData
    rect.set_transform(transform_)
    ax.add_patch(rect)
    return rect


def get_endpoint(rmodel, rdata, q: np.ndarray, tool_id: int):
    pin.framesForwardKinematics(rmodel, rdata, q)
    return rdata.oMf[tool_id].translation.copy()


def get_endpoint_traj(rmodel, rdata, xs: List[np.ndarray], tool_id: int):
    pts = []
    for i in range(len(xs)):
        pts.append(get_endpoint(rmodel, rdata, xs[i][: rmodel.nq], tool_id))
    return np.array(pts)


def compute_quasistatic(model: pin.Model, data: pin.Data, x0, acc):
    nq = model.nq
    q0 = x0[:nq]
    v0 = x0[nq:]
    return pin.rnea(model, data, q0, v0, acc)


def create_cartpole(N):
    model = pin.Model()
    geom_model = pin.GeometryModel()

    parent_id = 0

    cart_radius = 0.1
    cart_length = 5 * cart_radius
    cart_mass = 1.0
    joint_name = "joint_cart"

    geometry_placement = pin.SE3.Identity()
    geometry_placement.rotation = pin.Quaternion(
        np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])
    ).toRotationMatrix()

    joint_id = model.addJoint(
        parent_id, pin.JointModelPY(), pin.SE3.Identity(), joint_name
    )

    body_inertia = pin.Inertia.FromCylinder(cart_mass, cart_radius, cart_length)
    body_placement = geometry_placement
    model.appendBodyToJoint(
        joint_id, body_inertia, body_placement
    )  # We need to rotate the inertia as it is expressed in the LOCAL frame of the geometry

    shape_cart = fcl.Cylinder(cart_radius, cart_length)

    geom_cart = pin.GeometryObject(
        "shape_cart", joint_id, geometry_placement, shape_cart
    )
    geom_cart.meshColor = np.array([1.0, 0.1, 0.1, 1.0])
    geom_model.addGeometryObject(geom_cart)

    parent_id = joint_id
    joint_placement = pin.SE3.Identity()
    body_mass = 0.1
    body_radius = 0.1
    for k in range(N):
        joint_name = "joint_" + str(k + 1)
        joint_id = model.addJoint(
            parent_id, pin.JointModelRX(), joint_placement, joint_name
        )

        body_inertia = pin.Inertia.FromSphere(body_mass, body_radius)
        body_placement = joint_placement.copy()
        body_placement.translation[2] = 1.0
        model.appendBodyToJoint(joint_id, body_inertia, body_placement)

        geom1_name = "ball_" + str(k + 1)
        shape1 = fcl.Sphere(body_radius)
        geom1_obj = pin.GeometryObject(geom1_name, joint_id, body_placement, shape1)
        geom1_obj.meshColor = np.ones((4))
        geom_model.addGeometryObject(geom1_obj)

        geom2_name = "bar_" + str(k + 1)
        shape2 = fcl.Cylinder(body_radius / 4.0, body_placement.translation[2])
        shape2_placement = body_placement.copy()
        shape2_placement.translation[2] /= 2.0

        geom2_obj = pin.GeometryObject(geom2_name, joint_id, shape2_placement, shape2)
        geom2_obj.meshColor = np.array([0.0, 0.0, 0.0, 1.0])
        geom_model.addGeometryObject(geom2_obj)

        parent_id = joint_id
        joint_placement = body_placement.copy()
    end_frame = pin.Frame(
        "end_effector_frame",
        model.getJointId("joint_" + str(N)),
        0,
        body_placement,
        pin.FrameType(3),
    )
    model.addFrame(end_frame)
    geom_model.collision_pairs = []
    model.qinit = np.zeros(model.nq)
    model.qinit[1] = 0.0 * np.pi
    model.qref = pin.neutral(model)
    data = model.createData()
    geom_data = geom_model.createData()
    ddl = np.array([0])
    return model, geom_model, data, geom_data, ddl


def make_npendulum(N, ub=True, lengths=None):
    model = pin.Model()
    geom_model = pin.GeometryModel()

    parent_id = 0

    base_radius = 0.08
    shape_base = fcl.Sphere(base_radius)
    geom_base = pin.GeometryObject("base", 0, shape_base, pin.SE3.Identity())
    geom_base.meshColor = np.array([1.0, 0.1, 0.1, 1.0])
    geom_model.addGeometryObject(geom_base)

    joint_placement = pin.SE3.Identity()
    body_mass = 1.0
    body_radius = 0.06
    if lengths is None:
        lengths = [1.0 for _ in range(N)]

    for k in range(N):
        joint_name = "joint_" + str(k + 1)
        if ub:
            jmodel = pin.JointModelRUBX()
        else:
            jmodel = pin.JointModelRX()
        joint_id = model.addJoint(parent_id, jmodel, joint_placement, joint_name)

        body_inertia = pin.Inertia.FromSphere(body_mass, body_radius)
        body_placement = joint_placement.copy()
        body_placement.translation[2] = lengths[k]
        model.appendBodyToJoint(joint_id, body_inertia, body_placement)

        geom1_name = "ball_" + str(k + 1)
        shape1 = fcl.Sphere(body_radius)
        geom1_obj = pin.GeometryObject(geom1_name, joint_id, shape1, body_placement)
        geom1_obj.meshColor = np.ones((4))
        geom_model.addGeometryObject(geom1_obj)

        geom2_name = "bar_" + str(k + 1)
        shape2 = fcl.Cylinder(body_radius / 4, body_placement.translation[2])
        shape2_placement = body_placement.copy()
        shape2_placement.translation[2] /= 2.0

        geom2_obj = pin.GeometryObject(geom2_name, joint_id, shape2, shape2_placement)
        geom2_obj.meshColor = np.array([0.0, 0.0, 0.0, 1.0])
        geom_model.addGeometryObject(geom2_obj)

        parent_id = joint_id
        joint_placement = body_placement.copy()

    return model, geom_model, geom_model


def load_talos_upper_body():
    robot = erd.load("talos")
    qref = robot.model.referenceConfigurations["half_sitting"]
    locked_joints = list(range(1, 14))
    locked_joints += [23, 31]
    locked_joints += [32, 33]
    red_bot = robot.buildReducedRobot(locked_joints, qref)
    return red_bot


def add_namespace_prefix_to_models(model, collision_model, visual_model, namespace):
    """
    Lifted from this GitHub discussion:
    https://github.com/stack-of-tasks/pinocchio/discussions/1841
    """
    # Rename geometry objects in collision model:
    for geom in collision_model.geometryObjects:
        geom.name = f"{namespace}/{geom.name}"

    # Rename geometry objects in visual model:
    for geom in visual_model.geometryObjects:
        geom.name = f"{namespace}/{geom.name}"

    # Rename frames in model:
    for f in model.frames:
        f.name = f"{namespace}/{f.name}"

    # Rename joints in model:
    for k in range(len(model.names)):
        model.names[k] = f"{namespace}/{model.names[k]}"


def plot_controls_traj(
    times,
    us,
    ncols=2,
    axes=None,
    effort_limit=None,
    joint_names=None,
    rmodel: pin.Model = None,
):
    t0 = times[0]
    tf = times[-1]
    us = np.asarray(us)
    nu = us.shape[1]
    nrows, r = divmod(nu, ncols)
    nrows += int(r > 0)
    if axes is None:
        fig, axes = plt.subplots(nrows, ncols, sharex="col", figsize=(6.4, 6.4))
    else:
        fig = axes[0].get_figure()

    if rmodel is not None:
        effort_limit = rmodel.effortLimit
        joint_names = rmodel.names

    axes = axes.flatten()
    for i in range(nu):
        ax: plt.Axes = axes[i]
        ax.step(times[:-1], us[:, i])
        if effort_limit is not None:
            ylim = ax.get_ylim()
            ax.hlines(-effort_limit[i], t0, tf, colors="k", linestyles="--")
            ax.hlines(+effort_limit[i], t0, tf, colors="r", linestyles="dashdot")
            ax.set_ylim(*ylim)
        if joint_names is not None:
            joint_name = joint_names[i].lower()
            ax.set_ylabel(joint_name)
    fig.supxlabel("Time $t$")
    fig.suptitle("Controls trajectory")
    fig.tight_layout()
    return fig


def plot_velocity_traj(times, vs, rmodel: pin.Model, ncols=2):
    vs = np.asarray(vs)
    nv = vs.shape[1]
    idx_to_joint_id_map = {}
    jid = 0
    for i in range(nv):
        if i in rmodel.idx_vs.tolist():
            jid += 1
        idx_to_joint_id_map[i] = jid
    print(idx_to_joint_id_map)
    nrows, r = divmod(nv, ncols)
    nrows += int(r > 0)

    vel_limit = rmodel.velocityLimit
    t0 = times[0]
    tf = times[-1]

    fig, axes = plt.subplots(nrows, ncols)
    fig: plt.Figure
    axes = axes.flatten()
    for i in range(nv):
        ax: plt.Axes = axes[i]
        ax.plot(times, vs[:, i])
        jid = idx_to_joint_id_map[i]
        joint_name = rmodel.names[jid].lower()
        ylim = ax.get_ylim()
        ax.hlines(-vel_limit[i], t0, tf, colors="k", linestyles="--")
        ax.hlines(+vel_limit[i], t0, tf, colors="r", linestyles="dashdot")
        ax.set_ylim(*ylim)
        ax.set_ylabel(joint_name)

    fig.supxlabel("Time $t$")
    fig.suptitle("Velocity trajectory")
    fig.tight_layout()
    return fig


def underactuated_inv_dyn(
    model: pin.Model,
    data: pin.Data,
    q,
    v,
    B: np.ndarray,
    rcm: pin.RigidConstraintModel,
):
    nu = B.shape[1]

    pin.computeAllTerms(model, data, q, v)
    nle = data.nle.copy()
    print("rhs={}".format(nle))

    # 6x12
    rcd = rcm.createData()
    J = pin.getConstraintJacobian(model, data, rcm, rcd)

    mat = np.hstack([B, -J.T])
    ret = np.linalg.lstsq(mat, nle, rcond=None)
    err = mat @ ret[0] - nle
    err_norm = np.linalg.norm(err, np.inf)
    print("err norm = {}".format(err_norm))

    tau, fc = np.split(ret[0], [nu])
    return tau, fc
