#!/usr/bin/env python

from setuptools import setup, find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="uamc-spectrophore",
    version=get_version("spectrophore/__init__.py"),
    entry_points={
        "console_scripts": [
            "spectrophore=spectrophore.__main__:main",
        ],
    },
    author="Hans De Winter",
    author_email="hans.dewinter@uantwerpen.be",
    description=("Python implementation of the spectrophore descriptor"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/UAMCAntwerpen/spectrophore",
    download_url="https://github.com/silicos-it/spectrophore/releases/latest",
    packages=["spectrophore"],
    package_dir={"spectrophore": "spectrophore"},
    keywords=["uamc", "spectrophore", "rdkit", "cheminformatics"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: Freeware",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    install_requires=[
        "numba",
        "rdkit",
        "scipy",
        "numpy",
        "tqdm",
    ],
    extras_require={
        "develop": [
            "pytest",
            "pytest-cov",
            "ruff",
        ],
    },
    python_requires=">=3.9",
)
