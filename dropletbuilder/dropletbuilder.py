import mbuild
import numpy as np


def get_height(r, theta):
    """
    Helper function to get the height of a spherical cap
    """
    return r - r * np.cos(theta * np.pi / 180)


class Droplet(mbuild.Compound):
    """
    Builds a droplet on a lattice.

    Parameters
    ----------
    radius : int, default = 2
        radius of the droplet in nm
    angle : float, default = 90.0
        contact angle of the droplet in degrees
    fluid : mbuild.Compound or list of mbuild.Compound
        compounds to fill the droplet with
    density: float or list of float
        target density for the droplet in kg/m^3
    lattice: mbuild.Lattice
        lattice to build droplet on
    lattice_compound: mbuild.Compound
        compound to build lattice with
    x : float
        dimension of graphene sheet in x direction in nm
    y : float
        dimension of graphene sheet in y direction in nm
    NOTE: length of `fluid` must match length of `density`

    Attributes
    ----------
    see mbuild.Compound

    """

    def __init__(self, radius=2, angle=90.0, fluid=None, density=None,
                lattice=None, lattice_compound=None, x=None, y=None):

        super(Droplet, self).__init__()

        if fluid is None:
            raise ValueError('Fluid droplet compounds must be specified')
        if density is None:
            raise ValueError('Fluid density must be specified (units kg/m^3)')

        if x:
            if x < radius * 4:
                raise ValueError(
                    'Dimension x of sheet must be at least radius * 4')
            elif x > 100:
                raise ValueError(
                    'Dimension x of sheet must be less than 100 nm')
        else:
            x = radius * 4

        if y:
            if y < radius * 4:
                raise ValueError(
                    'Dimension y of sheet must be at least radius * 4')
            elif y > 100:
                raise ValueError(
                    'Dimension y of sheet must be less than 100 nm')
        else:
            y = radius * 4

        # Default to graphene lattice
        if lattice is None:

            if lattice_compound is not None:
                raise ValueError(
                    'If Lattice is None, defaults to a Graphene surface. ' +
                    'In this case, do not specify lattice_compound.'
                )

            lattice_compound = mbuild.Compound(name='C')
            lattice_spacing = [0.2456, 0.2456, 0.335]
            angles = [90.0, 90.0, 120.0]
            carbon_locations = [[0, 0, 0], [2 / 3, 1 / 3, 0]]
            basis = {lattice_compound.name: carbon_locations}
            lattice = mbuild.Lattice(
                lattice_spacing=lattice_spacing,
                angles=angles,
                lattice_points=basis)
            compound_dict = {lattice_compound.name: lattice_compound}

            factor = np.cos(np.pi / 6) # fixes non-cubic lattice
            # Estimate the number of lattice repeat units
            replicate = [int(x / 0.2456), int(y / 0.2456) * (1 / factor)]

            lat = lattice.populate(
                compound_dict=compound_dict,
                x=replicate[0],
                y=replicate[1],
                z=3
            )

            for particle in lat.particles():
                if particle.xyz[0][0] < 0:
                    particle.xyz[0][0] += lat.periodicity[0]
            lat.periodicity[1] *= factor

        else:
            if lattice_compound is None:
                raise ValueError('Lattice compounds must be specified')

            if not np.all(lattice.angles == 90.0):
                raise ValueError(
                    'Currently, only cubic lattices are supported. ' +
                    'If using Graphene, do not pass in a Lattice.'
                )

            compound_dict = {lattice_compound.name: lattice_compound}
            lat = lattice.populate(
                compound_dict=compound_dict,
                x=int(x/lattice.lattice_spacing[0]),
                y=int(y/lattice.lattice_spacing[1]),
                z=int(1.5/lattice.lattice_spacing[2]))

        sheet = mbuild.clone(lat)
        self.surface_height = np.max(sheet.xyz, axis=0)[2]
        coords = list(sheet.periodicity)

        height = get_height(radius, angle)
        sphere_coords = [coords[0] / 2, coords[1] / 2, radius, radius]
        sphere = mbuild.fill_sphere(
            compound=fluid, sphere=sphere_coords, density=density)

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

        sheet.name = 'LAT'
        sphere.name = 'FLD'
        sphere.xyz -= [0, 0, np.min(sphere.xyz, axis=0)[2]]
        sphere.xyz += [0, 0, self.surface_height + 0.3]

        self.add(sheet)
        self.add(sphere)
        self.periodicity[0] = sheet.periodicity[0]
        self.periodicity[1] = sheet.periodicity[1]
        self.periodicity[2] = radius * 5
