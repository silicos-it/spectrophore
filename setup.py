#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
	long_description = fh.read()

setup(
    name = "uamc-spectrophore",
    version = "1.0.1",
    author = "Hans De Winter",
    author_email = "hans.dewinter@uantwerpen.be",
    description = ("Python implementation of the spectrophore descriptor"),
    long_description = long_description,
	long_description_content_type="text/markdown",
    url = "https://github.com/UAMCAntwerpen/spectrophore",
	download_url = "https://github.com/UAMCAntwerpen/spectrophore/releases/tag/1.0.1",
    packages = find_packages(include=['spectrophore']),
	keywords = ['uamc', 'spectrophore', 'rdkit', 'cheminformatics'],
    classifiers = [
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
    install_requires = [
        "scipy",
        "numpy",
        "progressbar",
    ],
	python_requires = ">=3.6",
)
