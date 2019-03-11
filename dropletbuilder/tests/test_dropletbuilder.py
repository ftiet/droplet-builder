"""
Unit and regression test for the dropletbuilder package.
"""

# Import package, test suite, and other packages as needed
import dropletbuilder
import pytest
import sys

def test_dropletbuilder_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dropletbuilder" in sys.modules
