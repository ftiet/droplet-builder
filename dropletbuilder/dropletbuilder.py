import os
import mbuild as mb
import numpy as np
from foyer import Forcefield 
from pkg_resources import resource_filename

"""
dropletbuilder.py
Droplet on graphene builder

Handles the primary functions
"""

def get_fn(name):
    """
    Get the full path to one of the reference files shipped for utils.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the utils/ folder).

    """
    fn = resource_filename('dropletbuilder', os.path.join('utils', name))
    if not os.path.exists(fn):
        raise IOError('{} does not exist.'.format(fn))
    return fn

def get_height(r, theta):
    """
    Helper function to get the height of a spherical cap
    """
    return r - r*np.cos(theta*np.pi/180)

class GrapheneDroplet(mb.Compound):
    """
    Builds a droplet on a ~30x30 nm graphene sheet.

    Parameters
    ----------
    radius : int, default = 5
        radius of the droplet in nm
    angle : float, default = 90.0
        contact angle of the droplet in degrees

    Attributes
    ----------
    see mbuild.Compound

    """

    def __init__(self, radius=5, angle=90.0):
        super(GrapheneDroplet, self).__init__()

        water = mb.load(get_fn('tip3p.mol2'))
        graphene = mb.load(get_fn('graphene_sheet.pdb'))
        coords = list(graphene.periodicity)

        height = get_height(radius, angle)
        self.surface_height = np.max(graphene.xyz, axis=0)[2]
        sphere = mb.fill_sphere(compound=[water], sphere=[coords[0], coords[1], radius, radius], n_compounds=[2000])

        to_remove = []
        for child in sphere.children:
            for atom_coords in child.xyz:
                if height > radius:
                    if atom_coords[2] < height - radius:
                        to_remove += child
                        break
                else:
                    if atom_coords[2] < height:
                        to_remove += child
                        break

        sphere.remove(to_remove)

        sphere.name = 'H2O'
        sphere.xyz -= [0, 0, np.min(sphere.xyz, axis=0)[2]]
        sphere.xyz += [0, 0, self.surface_height + 0.3]

        sheet1 = mb.clone(graphene)
        sheet2 = mb.clone(graphene)
        sheet3 = mb.clone(graphene)
        sheet4 = mb.clone(graphene)
        sheet2.translate([coords[0],0,0])
        sheet3.translate([0,coords[1],0])
        sheet4.translate([coords[0],coords[1],0])

        sheet1.name = 'GPH'
        sheet2.name = 'GPH'
        sheet3.name = 'GPH'
        sheet4.name = 'GPH'

        self.add(sphere)
        self.add(sheet1)
        self.add(sheet2)
        self.add(sheet3)
        self.add(sheet4)
        self.periodicity[0] = coords[0]*2
        self.periodicity[1] = coords[1]*2
        self.periodicity[2] = 20
