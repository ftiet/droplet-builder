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
    def GrapheneDroplet(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        return GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997)

    @pytest.fixture
    def GrapheneDropletWithDims(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        return GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=4)


"""
Unit Tests for GrapheneDroplet class.
"""

class TestDropletBuilder(BaseTest):
    def test_dropletbuilder_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "dropletbuilder" in sys.modules

    def test_init_with_missing_fluid(self):
        with pytest.raises(ValueError, match="Fluid droplet compounds"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            GrapheneDroplet(radius=1, angle=90.0, density=997)
    
    def test_init_with_missing_density(self):
        with pytest.raises(ValueError, match="Fluid density"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            water = mbuild.load(get_fn('tip3p.mol2'))
            GrapheneDroplet(radius=1, angle=90.0, fluid=water)

    def test_init_with_too_small_x(self):
        with pytest.raises(ValueError, match="x .* at least"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            water = mbuild.load(get_fn('tip3p.mol2'))
            GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=1, y=4)

    def test_init_with_too_small_y(self):
        with pytest.raises(ValueError, match="y .* at least"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            water = mbuild.load(get_fn('tip3p.mol2'))
            GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=1)
    
    def test_init_with_too_large_x(self):
        with pytest.raises(ValueError, match="x .* 100"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            water = mbuild.load(get_fn('tip3p.mol2'))
            GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=101, y=4)

    def test_init_with_too_large_y(self):
        with pytest.raises(ValueError, match="y .* 100"):
            from dropletbuilder.dropletbuilder import GrapheneDroplet
            water = mbuild.load(get_fn('tip3p.mol2'))
            GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=4, y=101)

    def test_save(self, GrapheneDroplet):
        GrapheneDroplet.save('droplet.gro', overwrite=True, combine='all')

    def test_save_with_dims(self, GrapheneDropletWithDims):
        GrapheneDropletWithDims.save('droplet-with-dims.gro', overwrite=True, combine='all')

    def test_hierarchy(self, GrapheneDroplet):
        assert len(GrapheneDroplet.children) == 2

    def test_hierarchy_with_dims(self, GrapheneDropletWithDims):
        assert len(GrapheneDropletWithDims.children) == 2
    
    def test_lateral_dims_near_x_spec(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        x = 4
        droplet = GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=x)
        assert abs(droplet.boundingbox.lengths[0] - x) < 0.5

    def test_lateral_dims_near_y_spec(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        y = 4
        droplet = GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, y=y)
        assert abs(droplet.boundingbox.lengths[1] - y) < 0.5

    def test_lateral_dims_near_x_y_spec(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        water = mbuild.load(get_fn('tip3p.mol2'))
        x = 4
        y = 4
        droplet = GrapheneDroplet(radius=1, angle=90.0, fluid=water, density=997, x=x, y=y)
        assert (abs(droplet.boundingbox.lengths[0] - x) < 0.5 
                    and abs(droplet.boundingbox.lengths[1] - y) < 0.5)

    def test_lateral_dims_in_box(self, GrapheneDroplet):
        for child in GrapheneDroplet.children:
            if (np.min(child.xyz, axis=0)[0] < 0
                    or np.min(child.xyz, axis=0)[1] < 0):
                assert False
            if (np.max(child.xyz, axis=0)[0] > GrapheneDroplet.periodicity[0] or
                    np.max(child.xyz, axis=0)[1] > GrapheneDroplet.periodicity[1]):
                assert False
        assert True

    def test_lateral_dims_with_x_y_in_box(self, GrapheneDropletWithDims):
        for child in GrapheneDropletWithDims.children:
            if (np.min(child.xyz, axis=0)[0] < 0
                    or np.min(child.xyz, axis=0)[1] < 0):
                assert False
            if ((np.max(child.xyz, axis=0)[0] >
                 GrapheneDropletWithDims.periodicity[0])
                    or (np.max(child.xyz, axis=0)[1] >
                        GrapheneDropletWithDims.periodicity[1])):
                assert False
        assert True

    def test_n_fluid_particles(self, GrapheneDroplet):
        n_fluid_particles = 0
        for child in GrapheneDroplet.children:
            if child.name != 'LAT':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_n_fluid_particles_with_x_y(self, GrapheneDropletWithDims):
        n_fluid_particles = 0
        for child in GrapheneDropletWithDims.children:
            if child.name != 'LAT':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_fluid_particles_in_sheets(self, GrapheneDroplet):
        for child in GrapheneDroplet.children:
            if child.name != 'LAT':
                if (np.min(child.xyz, axis=0)[2] <
                        GrapheneDroplet.surface_height + 0.001):
                    assert False
        assert True

    def test_fluid_particles_in_sheets_with_x_y(self, GrapheneDropletWithDims):
        for child in GrapheneDropletWithDims.children:
            if child.name != 'LAT':
                if (np.min(child.xyz, axis=0)[2] <
                        GrapheneDropletWithDims.surface_height + 0.001):
                    assert False
        assert True
