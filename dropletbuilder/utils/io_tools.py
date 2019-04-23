import os
from pkg_resources import resource_filename


def get_fn(name):
    """
    Get the full path to one of the reference files shipped for utils.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the utils/ folder).
    """
    fn = resource_filename('dropletbuilder', os.path.join('utils/files', name))
    if not os.path.exists(fn):
        raise IOError('{} does not exist.'.format(fn))
    return fn