"""
Some useful utilities for arbow
"""

import subprocess


class IQtree:
    """
    Detect the iqtree version running and build appropriate
    command line for it.
    """
    def __init__(self, iqtree_path, iqtree_version):
        self.path = iqtree_path
        self.version = iqtree_version

    def _iqtree_version(self):
        pass


