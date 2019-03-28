import pytest
import sys
import numpy as np


class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def GrapheneDroplet(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        return GrapheneDroplet(radius=1, angle=90.0)

    @pytest.fixture
    def GrapheneDropletWithDims(self):
        from dropletbuilder.dropletbuilder import GrapheneDroplet
        return GrapheneDroplet(radius=1, angle=90.0, x=4, y=4)


"""
Unit Tests for GrapheneDroplet class.
"""

class TestDropletBuilder(BaseTest):
    def test_dropletbuilder_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "dropletbuilder" in sys.modules

    def test_save(self, GrapheneDroplet):
        GrapheneDroplet.save('droplet.gro', overwrite=True, combine='all')

    def test_save_with_dims(self, GrapheneDropletWithDims):
        GrapheneDropletWithDims.save('droplet-with-dims.gro', overwrite=True, combine='all')

    def test_hierarchy(self, GrapheneDroplet):
        assert len(GrapheneDroplet.children) == 2

    def test_hierarchy_with_dims(self, GrapheneDropletWithDims):
        assert len(GrapheneDropletWithDims.children) == 2

    def test_lateral_dims(self, GrapheneDroplet):
        for child in GrapheneDroplet.children:
            if (np.min(child.xyz, axis=0)[0] < 0
                    or np.min(child.xyz, axis=0)[1] < 0):
                assert False
            if (np.max(child.xyz, axis=0)[0] > GrapheneDroplet.periodicity[0] or
                    np.max(child.xyz, axis=0)[1] > GrapheneDroplet.periodicity[1]):
                assert False
        assert True

    def test_lateral_dims_with_dims(self, GrapheneDropletWithDims):
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
            if child.name != 'GPH':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_n_fluid_particles_with_dims(self, GrapheneDropletWithDims):
        n_fluid_particles = 0
        for child in GrapheneDropletWithDims.children:
            if child.name != 'GPH':
                n_fluid_particles += child.n_particles
        assert n_fluid_particles > 20 and n_fluid_particles < 150

    def test_fluid_particles_in_sheets(self, GrapheneDroplet):
        for child in GrapheneDroplet.children:
            if child.name != 'GPH':
                if (np.min(child.xyz, axis=0)[2] <
                        GrapheneDroplet.surface_height + 0.001):
                    assert False
        assert True

    def test_fluid_particles_in_sheets_with_dims(self, GrapheneDropletWithDims):
        for child in GrapheneDropletWithDims.children:
            if child.name != 'GPH':
                if (np.min(child.xyz, axis=0)[2] <
                        GrapheneDropletWithDims.surface_height + 0.001):
                    assert False
        assert True
