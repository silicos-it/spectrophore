#!/usr/bin/python
 
from distutils.core import setup
from distutils.extension import Extension

from os.path import join as pjoin

# Where to find extensions
SPEC = 'spectrophores'

extensions = []
extensions.append(
    Extension(pjoin(SPEC, 'CRotate'),
              language = 'c++',
              sources = [pjoin(SPEC, 'CRotate_boost.cpp')],
              libraries = ["boost_python"]))

setup(name = "spectrophore",
      version = '1.0.1',
      description = 'Spectrophore class to be used as Python library',
      platforms = ['Linux', 'Unix'],
      author = 'Fabio Mendes, Hans de Winter',
      author_email = 'fabiomendes.farm@gmail.com, hans.dewinter@uantwerpen.be',
      url = 'https://www.uantwerpen.be/nl/onderzoeksgroep/uamc/',
      download_url = 'https://github.com/UAMCAntwerpen/spectrophores/archive/spectrophore-1.0.1.tar.gz',
      ext_modules = extensions, 
      py_modules = [pjoin(SPEC, 'spectrophore')],
      package_dir = {'spectrophores': 'spectrophores'},
      packages = ['spectrophores'],
      scripts = ['sdf2spectrophore.py'],
      classifiers = [
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: End Users/Desktop',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics'
      ]
)


