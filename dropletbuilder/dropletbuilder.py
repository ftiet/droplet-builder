import mbuild as mb
import numpy as np
from foyer import Forcefield 
from ilforcefields.utils.utils import get_il
from ilforcefields.utils.utils import get_ff

"""
dropletbuilder.py
Droplet on graphene builder

Handles the primary functions
"""

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

        cat = get_il('emim')
        an = get_il('tf2n')
        cmpnds = [cat, an, cat, an]
        n_cmpnds = [150, 150, 150, 150]

        graphene = mb.load('./utils/graphene_sheet.pdb')
        coords = list(graphene.periodicity)

        height = get_height(radius, angle)
        sheet_height = 0.682
        sphere = mb.fill_sphere(compound=cmpnds, sphere=[coords[0], coords[1], radius, radius], n_compounds=n_cmpnds)

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

        sphere.xyz -= [0, 0, np.min(sphere.xyz, axis=0)[2]]
        sphere.xyz += [0, 0, sheet_height + 0.2]

        sheet1 = mb.clone(graphene)
        sheet2 = mb.clone(graphene)
        sheet3 = mb.clone(graphene)
        sheet4 = mb.clone(graphene)
        sheet2.translate([coords[0],0,0])
        sheet3.translate([0,coords[1],0])
        sheet4.translate([coords[0],coords[1],0])

        GRAPHENE = Forcefield('./utils/graphene.xml')
        LOPES = get_ff('lopes')

        system = GRAPHENE.apply(sheet1) + GRAPHENE.apply(sheet2) + GRAPHENE.apply(sheet3) + GRAPHENE.apply(sheet4)
        system = system + LOPES.apply(sphere.to_parmed(residues=['emim', 'tf2n']))

        system.box[0] = coords[0]*10*2
        system.box[1] = coords[1]*10*2
        system.box[2] = 200
        system.save('droplet.gro', overwrite=True)
        system.save('droplet.top', overwrite=True, combine='all')
