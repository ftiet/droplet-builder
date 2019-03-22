"""
Unit and regression test for the dropletbuilder package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
import numpy as np

@pytest.fixture
def GrapheneDroplet():
    from dropletbuilder.dropletbuilder import GrapheneDroplet
    return GrapheneDroplet()

"""
Unit Tests for GrapheneDroplet class.
"""

def test_dropletbuilder_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dropletbuilder" in sys.modules

def test_save(GrapheneDroplet):
    GrapheneDroplet.save('droplet.gro', overwrite=True)

def test_particles_in_sheets(GrapheneDroplet):
    for child in GrapheneDroplet.children:
        if child.name != 'GPH':
            if np.min(child.xyz, axis=0)[2] < GrapheneDroplet.surface_height + 0.001:
                assert False
    assert True
