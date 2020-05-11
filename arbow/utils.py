"""
Some useful utilities for arbow
"""
import os
import re
import subprocess
import logging
import sys
from shutil import which


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class IQTree:
    """
    Detect the iqtree version running and build appropriate
    command line for it.
    """

    def __init__(self, iqtree_path, iqtree_exec, iqtree_version):
        self.exec = iqtree_exec
        self.path = iqtree_path
        self.version = int(iqtree_version)
        if isinstance(self.path, str):
            self.path = [self.path]
        self.iq_path = None
        self.version_semver = None
        self.get_abs_exec()

    def _iqtree_version(self):
        version_str = re.compile('\\d\\.\\d\\.\\d')
        cmd = f"{self.iq_path} --version"
        proc = subprocess.run(cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              encoding='utf8')
        version_semver = version_str.findall(proc.stdout)
        if len(version_semver) == 0:
            logger.error(f"Was not able to parse version string from {self.iq_path}")
            sys.exit(1)
        if len(version_semver) > 1:
            logger.warning(f"Found too many version strings, taking the first")
        self.version_semver = version_semver[0]

    def get_abs_exec(self):
        iq_path = [path for path in self.path if which(self.exec, path=path)]
        if len(iq_path) == 0:
            logger.error(f"Could not find {self.exec} in {':'.join(self.path)}")
            sys.exit(1)
        if len(iq_path) > 0:
            logger.warning(f"Found multiple installs of iqtree in {':'.join(iq_path)} "
                           f"Assuming the first one is the right one!")
        self.iq_path = iq_path[0]
        self.iq_path = os.path.join(self.iq_path, self.exec)

    def check_version(self):
        self._iqtree_version()
        major, *_ = [int(version) for version in self.version_semver.split(".")]
        if major == self.version:
            logger.info(f"Expected IQTree version {self.version} and found {self.version_semver}")
            return
        if major < self.version:
            logger.error(f"Expected IQTree version {self.version} and found {self.version_semver}. Please upgrade!")
            sys.exit(1)
        if major > self.version:
            logger.error(f"Expected IQTree version {self.version} and found {self.version_semver}. Please downgrade!")
            sys.exit(1)
        logger.error("Something went wrong when testing for IQTree version. "
                      "Please file an issue with https://github.com/MDU-PHL/arbow/issues")
        sys.exit(1)

    def __str__(self):
        return self.iq_path
