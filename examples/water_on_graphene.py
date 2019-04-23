from dropletbuilder.dropletbuilder import GrapheneDroplet
import mbuild
from dropletbuilder.utils.io_tools import get_fn
from foyer import Forcefield

# water compound
water = mbuild.load(get_fn('tip3p.mol2'))

# build the system
system = GrapheneDroplet(radius=2, angle=110.0, fluid=water, density=997)

# get and apply forcefields
GPH = Forcefield(get_fn('graphene.xml'))
TIP3P = Forcefield(get_fn('tip3p.xml'))

for child in system.children:
    if child.name == 'LAT':
        lattice_pmd = GPH.apply(child.to_parmed(residues='LAT'))
    elif child.name == 'FLD':
        water_pmd = TIP3P.apply(child.to_parmed(residues='FLD'))
    else:
        continue

periodicity = system.periodicity
system = lattice_pmd + water_pmd 
system.box[2] = periodicity[2] * 100
system.save('water_graphene_drop.gro', overwrite=True, combine='all')
system.save('water_graphene_drop.top', overwrite=True, combine='all')