"""
This function comes from Kevin Dalton, and can be used to add the crystallographic 
cell axes for a given object to the bottom left of the screen -- VMD-style.
"""

from pymol import cmd
import numpy as np
from chempy import cpv


def cellbasis(angles, edges):
    """
    For the unit cell with given angles and edge lengths calculate the
    basis transformation (vectors) as a 4x4 matrix
    """
    rad = np.deg2rad(angles)
    basis = np.identity(4)
    basis[0][1] = np.cos(rad[2])
    basis[1][1] = np.sin(rad[2])
    basis[0][2] = np.cos(rad[1])
    basis[1][2] = (np.cos(rad[0]) - basis[0][1] * basis[0][2]) / basis[1][1]
    basis[2][2] = np.sqrt(1 - basis[0][2] ** 2 - basis[1][2] ** 2)
    edges.append(1.0)
    return basis * edges  # numpy.array multiplication!


class PutCenterCallback(object):
    prev_v = None

    def __init__(self, name, corner=0):
        self.name = name
        self.corner = corner
        self.cb_name = cmd.get_unused_name("_cb")

    def load(self):
        cmd.load_callback(self, self.cb_name)

    def __call__(self):
        if self.name not in cmd.get_names("objects"):
            import threading

            threading.Thread(None, cmd.delete, args=(self.cb_name,)).start()
            return

        v = cmd.get_view()
        if v == self.prev_v:
            return
        self.prev_v = v

        t = v[12:15]

        if self.corner:
            vp = cmd.get_viewport()
            R_mc = [v[0:3], v[3:6], v[6:9]]
            off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
            if self.corner in [2, 3]:
                off_c[0] *= -1
            if self.corner in [3, 4]:
                off_c[1] *= -1
            off_m = cpv.transform(R_mc, off_c)
            t = cpv.add(t, off_m)

        z = -v[11] / 30.0
        m = [z, 0, 0, 0, 0, z, 0, 0, 0, 0, z, 0, t[0] / z, t[1] / z, t[2] / z, 1]
        cmd.set_object_ttt(self.name, m)


def cell_axes(
    object=None,
    a_color="0xd95f02",
    b_color="0x1b9e77",
    c_color="0x7570b3",
    name="cell_axes",
):
    """
    DESCRIPTION
        Draw arrows corresponding to the crystallographic axes.
        The default color palette is colorblind friendly but close to the familiar red, green, and blue for the a, b, and, c axes respectively.
        (See https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3 for details)
    USAGE
        cell_axes [ object, [, a_color [, b_color [, c_color, [ name]]]]
    ARGUMENTS
        object = string: name of object to take cell definition from
        a_color = string: color of a-axis {default: '0xd95f02'}
        b_color = string: color of b-axis {default: '0x1b9e77'}
        c_color = string: color of c-axis {default: '0x7570b3'}
        name = string: name of the cgo object to create {default: cell_axes}
    """
    from pymol import cgo

    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]

    basis = cellbasis(cell_angles, cell_edges)
    print(basis)

    cmd.set("auto_zoom", 0)

    w = 0.06  # cylinder width
    l = 0.75  # cylinder length
    h = 0.25  # cone hight
    d = w * 1.618  # cone base diameter
    obj = []

    def add_axis_arrrow(obj, r, rgb, eps=1e-5):
        r = np.where(np.isclose(r, 0.0, atol=eps), 0.0, r)
        obj.extend(
            [cgo.CYLINDER, 0.0, 0.0, 0.0, l * r[0], l * r[1], l * r[2], w, *rgb, *rgb]
        )
        obj.extend(
            [
                cgo.CONE,
                l * r[0],
                l * r[1],
                l * r[2],
                (h + l) * r[0],
                (h + l) * r[1],
                (h + l) * r[2],
                d,
                0.0,
                *rgb,
                *rgb,
                1.0,
                1.0,
            ]
        )
        return obj

    A, B, C = basis[:3, :3].T
    r = A / np.linalg.norm(A)
    rgb = cmd.get_color_tuple(a_color)
    obj = add_axis_arrrow(obj, r, rgb)

    r = B / np.linalg.norm(B)
    rgb = cmd.get_color_tuple(b_color)
    obj = add_axis_arrrow(obj, r, rgb)

    r = C / np.linalg.norm(C)
    rgb = cmd.get_color_tuple(c_color)
    obj = add_axis_arrrow(obj, r, rgb)

    PutCenterCallback(name, 1).load()
    cmd.load_cgo(obj, name)


cmd.extend("cell_axes", cell_axes)
# tab-completion of arguments
cmd.auto_arg[0]["cell_axes"] = [cmd.object_sc, "object", ""]
