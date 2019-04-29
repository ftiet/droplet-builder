from dropletbuilder.dropletbuilder import Droplet
import mbuild
from dropletbuilder.utils.io_tools import get_fn
from foyer import Forcefield

# water compound
water = mbuild.load(get_fn('tip3p.mol2'))
water.name = 'H2O'

# build the system
system = Droplet(radius=2, angle=110.0, fluid=water, density=997)

# get and apply forcefields
GPH = Forcefield(get_fn('graphene.xml'))
TIP3P = Forcefield(get_fn('tip3p.xml'))

for child in system.children:
    if child.name == 'LAT':
        lattice_pmd = GPH.apply(child.to_parmed(residues='C'))
    elif child.name == 'FLD':
        water_pmd = TIP3P.apply(child.to_parmed(residues='H2O'))
    else:
        continue

periodicity = system.periodicity
system = lattice_pmd + water_pmd 
system.box[2] = periodicity[2] * 100
system.save('water_graphene_drop.gro', overwrite=True, combine='all')
system.save('water_graphene_drop.top', overwrite=True, combine='all')