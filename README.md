# Spectrophores

This module contains python code to calculate `spectrophores` from molecules. It is using the `RDKit` toolkit and it can also make use of available GPU's by means of the `numba` toolkit.

The technology and its applications have been described in [*Journal of Cheminformatics* (2018) **10**, 9](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0268-9). The paper is also included in this distribution.


## Installation

#### 1. Installation of RDKit and Numba

We recommend to install both `RDKit` and `Numba` using `Anaconda`. If `conda` is not yet available on your system, you should install this first following the instructions on the [Anaconda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation). [`Numba`]() is an open source JIT compiler that translates a subset of Python and NumPy code into fast machine code. [`RDKit`](https://www.rdkit.org) is open source cheminformatics software that provides the code to work with molecules.

The easiest way to have everything installed is using `conda`. First create a suitable environment in which you will install the `spectrophore` technology:

```console
> conda create --name spectrophore python=3
```

This will install a new `conda` environment with Python 3.9 installed in it. Now activate this environment:

```console
> conda activate spectrophore
```

`Numba` and `rdkit` can now be installed as follows (make sure you first activated the `spectrophore` environment):

```console
> conda install numba
> conda install cudatoolkit
> conda install -c conda-forge rdkit
```

You can test the `rdkit` installation by opening a `python` session from your command-line (assuming you are still in the activated `spectrophore` environment) and typing the following:

```python
>>> from rdkit import Chem
>>> mol = Chem.MolFromSmiles("C1CCCC1")
>>> print(mol.GetNumAtoms())
5
```

Similarly, you can test the `numba` installation with this small `python` snippet:

```python
>>> import numba
>>> numba.__version__
'0.51.2'
```

With `numba -s` command you can also check whether you have a CUDA device installed (check for the section `__CUDA Information__`).




  Step 1 - download the 'spectrophore-1.0.1.tar.gz' file from the 'dist/' folder

  Step 2 - untar the downloaded file and chdir into the created directory:

    > tar xvf spectrophore-1.0.1.tar.gz
    > cd spectrophore-1.0.1

  Step 3 - install:

    > sudo python setup.py install

  The installation process installs the following:

  - module 'spectrophore.py' is installed in the default third-party modules directory
    (for example '/usr/local/lib/python2.7/dist-packages/');

  - script 'sdf2spectrophore.py' is installed in the default third-party script directory
    (for example '/usr/local/bin/')


## Documentation

  All required documentation is contained in two ipython notebook files:
  - 'documentation/tutorial.ipynb'
  - 'documentation/manual.ipynb'


#### Reference and citation

If you use the `spectrophore` technology in your own research work, please cite as follows:

Gladysz, R.; Mendes Dos Santos, F.; Langenaeker, W.;  Thijs, G.; Augustyns, K.; De Winter, H. (2018) 'Spectrophores as one-dimensional descriptors calculated from three-dimensional atomic properties: applications ranging from scaffold hopping to multi-target virtual screening', *J. Cheminformatics* **10**, 9.
