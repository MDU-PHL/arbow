"""
Setup for Arbo
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
from arbow.version import __version__ as version

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="arbow",  # Required
    version=version,  # Required
    description="Cultivate your MSA to get better trees",  # Optional
    long_description=long_description,  # Optional
    long_description_content_type="text/markdown",  # Optional (see note above)
    url="https://github.com/MDPHL/arbow",  # Optional
    author="Microbiological Diagnostic Unit Public Health Laboratory",  # Optional
    author_email="andersgs@gmail.com",  # Optional
    classifiers=[  # Optional
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    keywords="bioinformatics phylogenetics multiple sequence alignment msa",  # Optional
    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    # package_dir={"": "arbow"},  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    # py_modules=["arbow"],
    #
    packages=find_packages(where="."),  # Required
    python_requires=">=3.6, <4",
    install_requires=["pandas", "click", "numpy", "biopython"],  # Optional
    # extras_require={"dev": ["check-manifest"], "test": ["coverage"]},  # Optional
    # package_data={"sample": ["package_data.dat"]},  # Optional
    # data_files=[("my_data", ["data/data_file"])],  # Optional
    entry_points={"console_scripts": ["arbow=arbow.__main__:main"]},  # Optional
    project_urls={  # Optional
        "Bug Reports": "https://github.com/MDUPHL/arbow/issues",
        "Source": "https://github.com/MDUPHL/arbow/",
    },
)
