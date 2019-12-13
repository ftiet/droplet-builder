from dropletbuilder.dropletbuilder import Droplet
import mbuild
from mbuild.examples.alkane.alkane import Alkane 
from dropletbuilder.utils.io_tools import get_fn
from foyer import Forcefield

# build the lattice
lattice_compound = mbuild.Compound(name='Au')
lattice_spacing = [0.40788, 0.40788, 0.40788]
lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
gold_locations = [[0., 0., 0.], [.5, .5, 0.], [.5, 0., .5], [0, .5, .5]]
basis = {lattice_compound.name: gold_locations}
gold_lattice = mbuild.Lattice(
    lattice_spacing=lattice_spacing,
    lattice_vectors=lattice_vector,
    lattice_points=basis)

# hexane compound
hexane = Alkane(n=6)
hexane.name = 'HEX'

# build the system
system = Droplet(radius=3, angle=90, fluid=hexane, density=655,
            lattice=gold_lattice, lattice_compound=lattice_compound)

# get and apply forcefields
AU = Forcefield(get_fn('heinz2008.xml'))
OPLSAA = Forcefield(name='oplsaa')

for child in system.children:
    if child.name == 'LAT':
        lattice_pmd = AU.apply(child.to_parmed(residues='Au'))
    elif child.name == 'FLD':
        hex_pmd = OPLSAA.apply(child.to_parmed(residues='HEX'))
    else:
        continue

periodicity = system.periodicity
system = lattice_pmd + hex_pmd
system.box[2] = periodicity[2] * 100
system.save('hexane_gold_drop.gro', overwrite=True, combine='all')
system.save('hexane_gold_drop.top', overwrite=True, combine='all')
