# Spectrophores

This module contains python code to calculate `spectrophores` from molecules. It is using the `RDKit` and `numba` toolkits.

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

#### 2. Installing the spectrophore code

With the `spectrophore` environment still active, you can now easily install the `spectrophore` module using:

```console
> pip install uamc-spectrophore
```

Check the installation by opening a `python` session and entering:

```python
>>> from spectrophore import spectrophore
>>> spectrophore.__version__
'1.0.1'
```


## Basic usage

Once you have installed all required tools and the `uamc-spectrophore` package, you are ready to use the tool. In its most simple form, `spectrophores` can be calculated as follows:

```python
>>> from spectrophore import spectrophore
>>> from rdkit import Chem
>>> from rdkit.Chem import AllChem
>>> mol = Chem.MolFromSmiles("c1ncncc1")
>>> mol = Chem.AddHs(mol)
>>> AllChem.EmbedMolecule(mol)
0
>>> calculator = spectrophore.SpectrophoreCalculator()
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> calculator.calculate(mol)
array([ 2.79612877e+00,  1.67804299e+00,  2.09592537e+00,  1.56546211e+00,
        3.52116032e+00,  3.46776209e+00,  2.62841820e+00,  2.60869245e+00,
        2.51220245e+00,  3.93695578e+00,  3.80399099e+00,  2.38644040e+00,
        4.03731656e+00,  1.29366764e+00,  5.21337207e-01,  1.23693624e-02,
       -1.46207663e+00,  5.63430303e+00, -5.07664172e-01, -5.76536023e-01,
        2.05383248e+00, -1.47519328e+00,  5.24238115e+00,  2.35758470e+00,
        9.59050950e+01,  2.72358109e+01,  9.49020528e+00, -2.42285499e+00,
       -3.54603042e+01,  1.17303747e+02, -1.52150804e+01, -1.72464368e+01,
        3.35391951e+01, -3.54603042e+01,  1.06633861e+02,  5.19670971e+01,
        1.83474518e+00,  6.12021101e-01,  1.28260191e-01, -1.10206015e-01,
       -7.59982097e-01,  2.61648795e+00, -3.23001423e-01, -3.69624105e-01,
        7.91595688e-01, -7.59982097e-01,  2.39011308e+00,  1.15975591e+00])
```

In the example shown, the first three lines import the required modules: module `spectrophore` for the calculation of spectrophores, module `Chem` to generate a RDKit molecule from a smiles string, and module `AllChem` to generate a 3D-conformation from the molecule. Next, a molecule is created from a smiles string (line 4), and a conformation is then generated at the 6'th line after adding hydrogen atoms on the 5'th line. Finally, on lines 7 and 8, a  `SpectrophoreCalculator` object is generated and this object is then used to calculate a `spectrophore` descriptor (line 8), which consists in its default form of 4 * 12 numbers.

> Note: a few words on the shape of a `spectrophore`:
>
> Each `spectrophore` consists of a set of floating point numbers, and this set is always a multiple of 4. The actual number count depends on how stereochemistry is treated in the calculation of spectrophores; this is controlled by the `stereo()` method:
- `stereo("none")`: the total number of numbers in a `spectrophore` is 48 (4 * 12; default),
- `stereo("unique")`: the total number of numbers in a `spectrophore` is 72 (4 * 18),
- `stereo("mirror")`: the total number of numbers in a `spectrophore` is 72 (4 * 18),
- `stereo("all")`: the total number of numbers in a `spectrophore` is 144 (4 * 36).
>
> Whatever the actual number of `spectrophore` points, these are always calculated according the same atomic properties. For example, consider a `spectrophore` of 4*n* points, then these points represent the following:
- **Points 1 to *n***: representing the interaction energies between the **atomic partial charges** and each of the *n* boxes;
- **Points *n*+1 to 2*n***: representing the interaction energies between the **atomic lipophilicities** and each of the *n* boxes;
- **Points 2*n*+1 to 3*n***: representing the interaction energies between the **atomic shape deviations** and each of the *n* boxes;
- **Points 3*n*+1 to 4*n***: representing the interaction energies between the **atomic electrophilicities** and each of the *n* boxes.
>
> Please have a look at the original publication form more information about the way these interaction energies are calculated, and what the `stereo()` method actually means.

If a molecule contains more than one 3D-conformation, then one may specify which conformation should be used for the calculation of `spectrophores`. As an example, consider the following code:

```python
>>> aspirin = Chem.MolFromSmiles("CC(Oc1ccccc1C(O)=O)=O")
>>> cids = AllChem.EmbedMultipleConfs(aspirin, numConfs=3)
>>> print(len(cids))
3
>>> for i in range(len(cids)): calculator.calculate(aspirin, i)
...
array([ 2.76333116e+00,  1.88272478e+00,  2.41594407e+00,  1.94137378e+00,
        4.24987797e+00,  3.31974773e+00,  2.69228201e+00,  2.63096102e+00,
        2.50042479e+00,  3.39582797e+00,  3.38116752e+00,  2.91764355e+00,
        1.84611805e+01,  3.72834082e+00,  4.30004028e+00,  3.98508974e+00,
       -5.72096326e-01,  1.06385182e+01,  3.90295711e-01,  1.00220548e+00,
        4.90389286e+00, -6.85082335e-01,  1.05008383e+01,  5.36045725e+00,
        2.53150467e+02,  8.69746911e+01,  4.62911497e+01,  4.05631990e+01,
       -1.85953640e+01,  2.49890738e+02,  1.83589314e+00, -1.13820244e+00,
        1.05522169e+02, -1.91212970e+01,  2.45247412e+02,  1.34650560e+02,
        7.69015379e+00,  2.66531515e+00,  1.41281704e+00,  1.23772227e+00,
       -5.72965434e-01,  7.56310105e+00,  4.82452251e-02, -4.02931606e-02,
        3.15069834e+00, -5.87971255e-01,  7.41740150e+00,  4.12861256e+00])
array([ 2.68627622e+00,  1.83634338e+00,  2.32891112e+00,  1.88955203e+00,
        4.15781984e+00,  3.26182940e+00,  2.62326816e+00,  2.58397849e+00,
        2.44327492e+00,  3.36734898e+00,  3.34772731e+00,  2.85710289e+00,
        1.77853239e+01,  3.67164582e+00,  4.15154181e+00,  3.83242016e+00,
       -5.67838738e-01,  1.05490695e+01,  3.79239626e-01,  9.99403083e-01,
        4.89491453e+00, -6.78547495e-01,  1.04110889e+01,  5.28024565e+00,
        2.38736866e+02,  8.40744401e+01,  4.28507882e+01,  3.71381291e+01,
       -1.84811217e+01,  2.45651168e+02,  1.39771377e+00, -1.54160197e+00,
        1.03879971e+02, -1.89940228e+01,  2.40980957e+02,  1.31583317e+02,
        7.24363945e+00,  2.58010744e+00,  1.30258311e+00,  1.12788219e+00,
       -5.70649605e-01,  7.45225558e+00,  3.40850493e-02, -5.45057478e-02,
        3.10737509e+00, -5.85184094e-01,  7.30536160e+00,  4.04229515e+00])
array([ 2.91469475e+00,  1.93172622e+00,  2.43339676e+00,  1.94337237e+00,
        4.34287437e+00,  3.42072192e+00,  2.74090844e+00,  2.70491554e+00,
        2.55278863e+00,  3.53028262e+00,  3.52315311e+00,  2.97409546e+00,
        1.84625592e+01,  3.73065078e+00,  4.32493133e+00,  4.00600516e+00,
       -5.65341277e-01,  1.06834328e+01,  4.06627497e-01,  1.05469332e+00,
        4.96638708e+00, -6.75702313e-01,  1.05473789e+01,  5.34733478e+00,
        2.46382770e+02,  8.56744422e+01,  4.42393367e+01,  3.85624404e+01,
       -1.85359724e+01,  2.47716551e+02,  1.51078638e+00, -1.53106741e+00,
        1.04196592e+02, -1.90381597e+01,  2.43050756e+02,  1.33308450e+02,
        7.47357302e+00,  2.62967431e+00,  1.34786441e+00,  1.17423161e+00,
       -5.71793341e-01,  7.51242068e+00,  3.92631649e-02, -5.22361423e-02,
        3.11897607e+00, -5.86099546e-01,  7.36595668e+00,  4.09366717e+00])
```

One can easily visualise `spectrophores` by plotting the actual values. For example, consider the following snippet:

```python
>>> import matplotlib.pyplot as plt
>>> mol = Chem.MolFromSmiles("CC(CCC1=CC=CC=C1Cl)N1CCOCC1")
>>> mol = Chem.AddHs(mol)
>>> cids = AllChem.EmbedMultipleConfs(mol, numConfs = 10)
>>> spectrophores = []
>>> for cid in cids: spectrophores.append(calculator.calculate(mol, cid))
...
>>> for i in range(3): plt.plot(range(1,49), spectrophores[i], label='Conf %d' % (i+1))
...
[<matplotlib.lines.Line2D object at 0x7faf42154be0>]
[<matplotlib.lines.Line2D object at 0x7faf420db880>]
[<matplotlib.lines.Line2D object at 0x7faf420df5e0>]
>>> plt.legend(loc='upper left')
>>> plt.suptitle("Spectrophores calculated for three molecular conformations")
>>> plt.savefig("exampleplot1.png")
```

which generates the following plot:

![Three conformations](exampleplot1.png)

Similarly, one can easily compare the `spectrophores` from two different molecules, and quantify the difference:

```python
>>> plt.close()
>>> mols = [Chem.MolFromSmiles("ClC(Br)(I)F"), Chem.MolFromSmiles("CC(CCC1=CC=CC=C1Cl)N1CCOCC1")]
>>> for i in range(2):
...    mols[i] = Chem.AddHs(mols[i])
...    AllChem.EmbedMolecule(mols[i])
...    spectrophores.append(calculator.calculate(mols[i]))
...
0
0
>>> for i in range(2): plt.plot(range(1,49), spectrophores[i], label='Molecule %d' % (i+1))
...
[<matplotlib.lines.Line2D object at 0x7faf420df520>]
[<matplotlib.lines.Line2D object at 0x7faf420df4c0>]
>>> plt.savefig("exampleplot2.png")
>>> from scipy.spatial import distance
>>> distance.euclidean(spectrophores[0],spectrophores[1])
2162.512188773051
```

![Two molecules](exampleplot2.png)

From the last example, it is clear that the actual `spectrophore` values may differ a lot depending on the type of molecule. Also, the absolute values are depending on the property type, with some properties leading to large values (e.g. shape deviation) and others very small. For this reason, a number of normalisation methods are provided as shown below.


## Methods

### `resolution()`

The `resolution()` method controls the smallest distance between the molecule and the surrounding box. By default this value is set to 3.0 A. The `resolution()` can be specified at the moment of class creation, or later on using the `resolution()` method:

```python
>>> mol = Chem.MolFromSmiles("ClC(Br)(I)F")
>>> AllChem.EmbedMolecule(mol)
0
>>> calculator = spectrophore.SpectrophoreCalculator()  # Default of 3.0
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.30117366411301427
>>> calculator = spectrophore.SpectrophoreCalculator(resolution = 3.0)
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.30117366411301427
>>> calculator = spectrophore.SpectrophoreCalculator(resolution = 5.0)
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.07899008784089182
>>> calculator.resolution(10.0)
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.009178129128602913
```

The larger the resolution value (e.g. 5.0 versus 3.0 A), the smaller the interaction energies and corresponding `spectrophore` values.

Calling the `resolution()` method without an argument returns the current resolution value:

```python
>>> calculator.resolution()
10.0
```


### `accuracy()`

The `accuracy()` method controls the angular stepsize by which the molecule is rotated within the cages. By default this value is set to 20°. This parameter can be modified either at class creation, or using the `accuracy()` method later on. The accuracy should be an integer fraction of 360, hence 360 modulus *accuracy* should be equal to 0. The smaller the accuracy value (meaning smaller angular stepsizes), the longer the computation time:

```python
>>> calculator = spectrophore.SpectrophoreCalculator(accuracy = 20.0) # Default
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.2947654850479195
>>> calculator = spectrophore.SpectrophoreCalculator()
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.2947654850479195
>>> calculator = spectrophore.SpectrophoreCalculator(accuracy = 2.0)
Probes initialised: 48 number of probes in total
Only using 12 probes
>>> spec = calculator.calculate(mol)
>>> print(spec[0])
0.30442130872737927
```

Calling the `accuracy()` method without an argument returns the current accuracy value:

```python
>>> calculator.accuracy()
2
```



### `normalization()`

With the `normalization()` method, one can specify the type of `spectrophore` normalization. There are four possibilities:
- `normalization("none")`: no normalization is applied and the `spectrophore` values are the raw calculated interaction energies (multiplied by -100),
- `normalization("mean")`: for each property, the average value is calculated and each of the individual `spectrophore` property value are reduced by these mean values. This centers the calculated values around 0,
- `normalization("std")`: for each property, the standard deviation is calculated and each of the individual `spectrophore` property value is divided by these standard deviations,
- `normalization("all")`: each spectrophore value is normalized by mean and standard deviation.

The default value is "none", but a setting of "all" is also a good choice:

```python
>>> calculator.normalization("none")
>>> spec = calculator.calculate(mol)
>>> print(spec[:12])
[0.293261   0.45461287 0.4156776  1.10858653 2.80451603 1.23886813
 1.15383657 1.08827509 1.10487316 2.76474151 1.62110108 1.22473783]
>>> calculator.normalization("mean")
>>> spec = calculator.calculate(mol)
print(spec[:12])
[-0.97949629 -0.81814442 -0.85707968 -0.16417076  1.53175875 -0.03388915
 -0.11892071 -0.18448219 -0.16788413  1.49198423  0.3483438  -0.04801945]
>>> sum(spec[:12])
-3.774758283725532e-15
>>> calculator.normalization("std")
>>> spec = calculator.calculate(mol)
>>> print(spec[:12])
[0.37955444 0.58838486 0.53799271 1.43479337 3.62975819 1.60341096
 1.49335846 1.4085052  1.42998733 3.57827982 2.09811777 1.58512275]
>>> calculator.normalization("all")
>>> spec = calculator.calculate(mol)
>>> print(spec[:12])
[-1.26771772 -1.05888729 -1.10927945 -0.21247878  1.98248603 -0.04386119
 -0.15391369 -0.23876695 -0.21728483  1.93100767  0.45084561 -0.0621494 ]
>>> sum(spec[:12])
-5.044575868140555e-15
```

Using a normalization over 'all' makes it more easier to compare `spectrophores`  between molecules:

```python
>>> mols = [Chem.MolFromSmiles("ClC(Br)(I)F"), Chem.MolFromSmiles("CC(CCC1=CC=CC=C1Cl)N1CCOCC1")]
>>> spectrophores = []
>>> for i in range(2):
...    mols[i] = Chem.AddHs(mols[i])
...    AllChem.EmbedMolecule(mols[i])
...    spectrophores.append(calculator.calculate(mols[i]))
...
0
0
>>> for i in range(2): plt.plot(range(1,49), spectrophores[i], label='Molecule %d' % (i+1))
...
[<matplotlib.lines.Line2D object at 0x7faf420df520>]
[<matplotlib.lines.Line2D object at 0x7faf420df4c0>]
>>> plt.savefig("exampleplot3.png")
>>> from scipy.spatial import distance
>>> distance.euclidean(spectrophores[0],spectrophores[1])
9.072753455280592
```

![Two molecules](exampleplot3.png)

The same holds true when comparing `spectrophores` from different conformations:

```python
>>> plt.close()
>>> spectrophores = []
>>> mol = Chem.MolFromSmiles("CC(CCC1=CC=CC=C1Cl)N1CCOCC1")
>>> mol = Chem.AddHs(mol)
>>> cids = AllChem.EmbedMultipleConfs(mol, numConfs = 10)
>>> calculator.normalization("all")
>>> for cid in cids: spectrophores.append(calculator.calculate(mol, cid))
...
>>> for i in range(3): plt.plot(range(1,49), spectrophores[i], label='Conf %d' % (i+1))
...
[<matplotlib.lines.Line2D object at 0x7faf42154be0>]
[<matplotlib.lines.Line2D object at 0x7faf420db880>]
[<matplotlib.lines.Line2D object at 0x7faf420df5e0>]
>>> plt.legend(loc='upper left')
>>> plt.suptitle("Spectrophores calculated for three molecular conformations")
>>> plt.savefig("exampleplot4.png")
```

![Three conformations](exampleplot4.png)


### `stereo()`

The `stereo()` method specifies the kind of cages to be used. The reason for this is that some of the cages that are used to calculate `spectrophores` have a stereospecific distribution of the interaction points:

![Stereo cages](exampleplot5.png)

There are four possibilities:
- `stereo("none")`: no stereospecificity (default). `Spectrophores` are generated using cages that are not stereospecific. For most applications, these `spectrophores` will suffice,
- `stereo("unique")`: unique stereospecificity. `Spectrophores` are generated using unique stereospecific cages,
- `stereo("mirror")`: mirror stereospecificity. Mirror stereospecific `spectrophores` are `spectrophores` resulting from the mirror enantiomeric form of the input molecules,
- `stereo("all")`: all cages are used. This results in the longest `spectrophores` and should only in specific cases be used.

The differences between the corresponding data points of unique and mirror stereospecific `spectrophores` are very small and require very long calculation times to obtain a sufficiently high quality level. This increased quality level is triggered by the `accuracy` setting and will result in calculation times being increased by at least a factor 100. As a consequence, it is recommended to apply this increased accuracy only in combination with a limited number of molecules, and when the small differences between the stereospecific `spectrophores` are really critical. However, for the vast majority of virtual screening applications, this increased accuracy is not required as long as it is not the intention to draw conclusions about differences in the underlying molecular stereoselectivity. Non-stereospecific `spectrophores` will therefore suffice for most applications.


## Interpreting `spectrophores`

A `spectrophore` is a vector of real number and has a certain length. The length depends on the used `stereo` method and the number of properties. The standard setting uses a set of non-stereospecific probes in combination with four properties:
- property 1: atomic partial charges
- property 2: atomic lipophilicities
- property 3: atomic shape deviations
- property 4: atomic electrophilicties

The combination of four properties and the set of non-stereospecific probes leads to a `spectrophore` vector length of 48. The use of other probes leads to other vector lengths, as summarised in this table:

| Stereospecificity | Number of probes | Number of properties | Length |
| ----------------- |:----------------:|:--------------------:|:------:|
| none              | 12               | 4                    | 48     |
| unique            | 18               | 4                    | 72     |
| mirror            | 18               | 4                    | 72     |
| all               | 36               | 4                    | 144    |

The general layout of a `spectrophore`, irrespective of its length, is always:

| Property 1         | Property 2         |         Property 3 |         Property 4 |
|:------------------:|:------------------:|:------------------:|:------------------:|
| probe 1..probe *n* | probe 1..probe *n* | probe 1..probe *n* | probe 1..probe *n* |

meaning that the first *n* values (with *n* being the number of probes) are calculated using property 1 (partial charges), then another *n* values (*n*+1 up to 2*n*) calculated using property 2 (lipophilicities), and so forth.



## Reference and citation

If you use the `spectrophore` technology in your own research work, please cite as follows:

Gladysz, R.; Mendes Dos Santos, F.; Langenaeker, W.;  Thijs, G.; Augustyns, K.; De Winter, H. (2018) 'Spectrophores as one-dimensional descriptors calculated from three-dimensional atomic properties: applications ranging from scaffold hopping to multi-target virtual screening', *J. Cheminformatics* **10**, 9.
