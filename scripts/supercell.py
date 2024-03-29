"""
(c) 2020 Jack Greisman
Based on script by Thomas Holder (2010)
"""

from pymol import cmd, cgo, xray
import numpy as np


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


def supercell(a=1, b=1, c=1, object=None, color="blue", name="supercell", withmates=1):
    """
    DESCRIPTION

    Draw a supercell, as requested by Nicolas Bock on the pymol-users
    mailing list (Subject: [PyMOL] feature request: supercell construction
    Date: 04/12/2010 10:12:17 PM (Mon, 12 Apr 2010 14:12:17 -0600))

    USAGE

    supercell a, b, c [, object [, color [, name [, withmates]]]]

    ARGUMENTS

    a, b, c = integer: repeat cell in x,y,z direction a,b,c times
    {default: 1,1,1}

    object = string: name of object to take cell definition from

    color = string: color of cell {default: blue}

    name = string: name of the cgo object to create {default: supercell}

    withmates = bool: also create symmetry mates in displayed cells
    {default: 1}

    SEE ALSO

    show cell

    """
    if object is None:
        object = cmd.get_object_list()[0]
    withmates = int(withmates)

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]

    basis = cellbasis(cell_angles, cell_edges)
    assert isinstance(basis, np.ndarray)

    ts = list()
    for i in range(int(a)):
        for j in range(int(b)):
            for k in range(int(c)):
                ts.append([i, j, k])

    obj = [
        cgo.BEGIN,
        cgo.LINES,
        cgo.COLOR,
    ]
    obj.extend(cmd.get_color_tuple(color))

    for t in ts:
        shift = basis[0:3, 0:3] * t
        shift = shift[:, 0] + shift[:, 1] + shift[:, 2]

        for i in range(3):
            vi = basis[0:3, i]
            vj = [
                np.array([0.0, 0.0, 0.0]),
                basis[0:3, (i + 1) % 3],
                basis[0:3, (i + 2) % 3],
                basis[0:3, (i + 1) % 3] + basis[0:3, (i + 2) % 3],
            ]
            for j in range(4):
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j]).tolist())
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j] + vi).tolist())

            if withmates:
                symexpcell("m%d%d%d_" % tuple(t), object, *t)

    obj.append(cgo.END)

    cmd.delete(name)
    cmd.load_cgo(obj, name)


def symexpcell(prefix="mate", object=None, a=0, b=0, c=0):
    """
    DESCRIPTION

    Creates all symmetry-related objects for the specified object that
    occur with their bounding box center within the unit cell.

    USAGE

    symexpcell prefix, object, [a, b, c]

    ARGUMENTS

    prefix = string: prefix of new objects

    object = string: object for which to create symmetry mates

    a, b, c = integer: create neighboring cell {default: 0,0,0}

    SEE ALSO

    symexp, http://www.pymolwiki.org/index.php/SuperSym
    """
    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
    spacegroup = sym[6]

    basis = cellbasis(cell_angles, cell_edges)
    basis = np.matrix(basis)

    extent = cmd.get_extent(object)
    center = sum(np.array(extent)) * 0.5
    center = np.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center

    extra_shift = [[float(i)] for i in (a, b, c)]

    i = 0
    matrices = xray.sg_sym_to_mat_list(spacegroup)
    for mat in matrices:
        i += 1

        mat = np.matrix(mat)
        shift = np.floor(mat * center_cell)
        mat[0:3, 3] -= shift[0:3, 0]
        mat[0:3, 3] += extra_shift

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = "%s%d" % (prefix, i)
        cmd.create(name, object)
        cmd.transform_object(name, mat_list, 0)
        cmd.color(i + 1, name)


cmd.extend("symexpcell", symexpcell)
cmd.extend("supercell", supercell)

# tab-completion of arguments
cmd.auto_arg[3]["supercell"] = [cmd.object_sc, "object", ""]
