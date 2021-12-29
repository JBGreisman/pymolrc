'''
publication.py
(c) 2020 Jack Greisman

Publication settings
'''

from pymol import cmd, util

def publication():
    """
    DESCRIPTION

    Settings to render publication quality figures

    USAGE

    publication
    """
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", True)
    cmd.set("mesh_radius", 0.01)
    cmd.set("antialias", 1)
    cmd.set("direct", 0.5)
    util.ray_shadows('none')
    return

cmd.extend('publication', publication)
