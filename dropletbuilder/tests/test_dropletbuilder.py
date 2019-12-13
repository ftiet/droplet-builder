import pytest
import sys
import numpy as np
import mbuild

from dropletbuilder.utils.io_tools import get_fn


class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def GoldLattice(self):
        lattice_compound = mbuild.Compound(name='Au')
        lattice_spacing = [0.40788, 0.40788, 0.40788]
        lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        gold_locations = [[0., 0., 0.], [.5, .5, 0.], [.5, 0., .5], [0, .5, .5]]
        basis = {lattice_compound.name: gold_locations}
        gold_lattice = mbuild.Lattice(
            lattice_spacing=lattice_spacing,
            lattice_vectors=lattice_vector,
            lattice_points=basis)
        return gold_lattice

    @pytest.fixture
    def Droplet(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        return Droplet(radius=1, angle=90.0, fluid=water, density=997)

    @pytest.fixture
    def DropletWithDims(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        return Droplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=4)


"""
Unit Tests for Droplet class.
"""

class TestDropletBuilder(BaseTest):
    def test_dropletbuilder_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "dropletbuilder" in sys.modules

    def test_init_with_missing_fluid(self):
        from dropletbuilder.dropletbuilder import Droplet
        with pytest.raises(ValueError, match="Fluid droplet compounds"):
            Droplet(radius=1, angle=90.0, density=997)
    
    def test_init_with_missing_density(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="Fluid density"):
            Droplet(radius=1, angle=90.0, fluid=water)

    def test_init_without_lattice_with_lattice_compound(self):
        lattice_compound = mbuild.Compound(name='Au')
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="do not specify lattice_compound"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=4, lattice_compound=lattice_compound)

    def test_init_with_lattice_without_lattice_compound(self, GoldLattice):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="Lattice compounds"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=4, lattice=GoldLattice)

    def test_init_with_too_small_x(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="x .* at least"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=1, y=4)

    def test_init_with_too_small_y(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="y .* at least"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=1)
    
    def test_init_with_too_large_x(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="x .* 100"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=101, y=4)

    def test_init_with_too_large_y(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        with pytest.raises(ValueError, match="y .* 100"):
            Droplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=101)
    
    def test_save(self, Droplet):
        Droplet.save('droplet.gro', overwrite=True, combine='all')

    def test_save_with_dims(self, DropletWithDims):
        DropletWithDims.save('droplet-with-dims.gro', overwrite=True, combine='all')

    def test_hierarchy(self, Droplet):
        assert len(Droplet.children) == 2

    def test_hierarchy_with_dims(self, DropletWithDims):
        assert len(DropletWithDims.children) == 2
    
    def test_lateral_dims_near_x_spec(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        x = 4
        droplet = Droplet(radius=1, angle=90.0, fluid=water, density=997, x=x)
        assert abs(droplet.boundingbox.lengths[0] - x) < 0.5

    def test_lateral_dims_near_y_spec(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        y = 4
        droplet = Droplet(radius=1, angle=90.0, fluid=water, density=997, y=y)
        assert abs(droplet.boundingbox.lengths[1] - y) < 0.5

    def test_lateral_dims_near_x_y_spec(self):
        from dropletbuilder.dropletbuilder import Droplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        x = 4
        y = 4
        droplet = Droplet(radius=1, angle=90.0, fluid=water, density=997, x=x, y=y)
        assert (abs(droplet.boundingbox.lengths[0] - x) < 0.5 
                    and abs(droplet.boundingbox.lengths[1] - y) < 0.5)

    def test_lateral_dims_in_box(self, Droplet):
        for child in Droplet.children:
            if (np.min(child.xyz, axis=0)[0] < 0
                    or np.min(child.xyz, axis=0)[1] < 0):
                assert False
            if (np.max(child.xyz, axis=0)[0] > Droplet.periodicity[0] or
                    np.max(child.xyz, axis=0)[1] > Droplet.periodicity[1]):
                assert False
        assert True

    def test_lateral_dims_with_x_y_in_box(self, DropletWithDims):
        for child in DropletWithDims.children:
            if (np.min(child.xyz, axis=0)[0] < 0
                    or np.min(child.xyz, axis=0)[1] < 0):
                assert False
            if ((np.max(child.xyz, axis=0)[0] >
                 DropletWithDims.periodicity[0])
                    or (np.max(child.xyz, axis=0)[1] >
                        DropletWithDims.periodicity[1])):
                assert False
        assert True

    def test_n_fluid_particles(self, Droplet):
        n_fluid_particles = 0
        for child in Droplet.children:
            if child.name != 'LAT':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_n_fluid_particles_with_x_y(self, DropletWithDims):
        n_fluid_particles = 0
        for child in DropletWithDims.children:
            if child.name != 'LAT':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_fluid_particles_in_sheets(self, Droplet):
        for child in Droplet.children:
            if child.name != 'LAT':
                if (np.min(child.xyz, axis=0)[2] <
                        Droplet.surface_height + 0.001):
                    assert False
        assert True

    def test_fluid_particles_in_sheets_with_x_y(self, DropletWithDims):
        for child in DropletWithDims.children:
            if child.name != 'LAT':
                if (np.min(child.xyz, axis=0)[2] <
                        DropletWithDims.surface_height + 0.001):
                    assert False
        assert True
