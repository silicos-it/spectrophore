import pickle
import pytest
import os
import subprocess

import numpy as np

from spectrophore import spectrophore
from rdkit import Chem
from rdkit.Chem import AllChem


TOLERANCE = 5e-4


# Load reference mol
@pytest.fixture
def ref_mol():
    with open(os.path.join(os.path.dirname(__file__), "ref_mol.pkl"), "rb") as f:
        return pickle.load(f)


# Load reference data
@pytest.fixture
def ref_data():
    with open(os.path.join(os.path.dirname(__file__), "ref_data.pkl"), "rb") as f:
        return pickle.load(f)


def test_spectrophore_creation_default():
    """
    Test that a spectrophore calculator can be created
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert isinstance(calculator, spectrophore.SpectrophoreCalculator)


def test_spectrophore_creation_normalization_1():
    """
    Test the normalization default setting
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert calculator.normalization() == "all"


def test_spectrophore_creation_normalization_2():
    """
    Test the normalization 'all' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(normalization="all")
    assert calculator.normalization() == "all"


def test_spectrophore_creation_normalization_3():
    """
    Test the normalization 'none' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(normalization="none")
    assert calculator.normalization() == "none"


def test_spectrophore_creation_normalization_4():
    """
    Test the normalization 'std' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(normalization="std")
    assert calculator.normalization() == "std"


def test_spectrophore_creation_normalization_5():
    """
    Test the normalization 'mean' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(normalization="mean")
    assert calculator.normalization() == "mean"


def test_spectrophore_creation_stereo_1():
    """
    Test the stereo default setting
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert calculator.stereo() == "none"


def test_spectrophore_creation_stereo_2():
    """
    Test the stereo 'none' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(stereo="none")
    assert calculator.stereo() == "none"


def test_spectrophore_creation_stereo_3():
    """
    Test the stereo 'unique' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(stereo="unique")
    assert calculator.stereo() == "unique"


def test_spectrophore_creation_stereo_4():
    """
    Test the stereo 'mirror' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(stereo="mirror")
    assert calculator.stereo() == "mirror"


def test_spectrophore_creation_stereo_5():
    """
    Test the stereo 'all' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(stereo="all")
    assert calculator.stereo() == "all"


def test_spectrophore_creation_resolution_1():
    """
    Test the resolution default setting
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert calculator.resolution() == 3.0


def test_spectrophore_creation_resolution_2():
    """
    Test the resolution setting
    """
    calculator = spectrophore.SpectrophoreCalculator(resolution=5.0)
    assert calculator.resolution() == 5.0


def test_spectrophore_creation_accuracy_1():
    """
    Test the accuracy default setting
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert calculator.accuracy() == 20


def test_spectrophore_creation_accuracy_2():
    """
    Test the accuracy setting
    """
    calculator = spectrophore.SpectrophoreCalculator(accuracy=30)
    assert calculator.accuracy() == 30


def test_spectrophore_creation_mode_1():
    """
    Test the mode default setting
    """
    calculator = spectrophore.SpectrophoreCalculator()
    assert calculator.mode() == "classic"


def test_spectrophore_creation_mode_2():
    """
    Test the mode 'classic' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(mode="classic")
    assert calculator.mode() == "classic"


def test_spectrophore_creation_mode_3():
    """
    Test the mode 'full' setting
    """
    calculator = spectrophore.SpectrophoreCalculator(mode="full")
    assert calculator.mode() == "full"


def test_spectrophore_valueerror_natoms():
    """
    Test the error raised when less than 3 atoms are present in the molecule
    """
    mol = Chem.MolFromSmiles("Br")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    calculator = spectrophore.SpectrophoreCalculator()

    with pytest.raises(
        ValueError, match=">=3 atoms are needed in molecule, only 2 given"
    ):
        calculator.calculate(mol)


def test_spectrophore_valueerror_conformers():
    """
    Test the error raised when no conformation is available
    """
    mol = Chem.MolFromSmiles("c1ccc(CCN)cc1")
    calculator = spectrophore.SpectrophoreCalculator()
    with pytest.raises(ValueError, match="Conformer ID 0 is not a valid"):
        calculator.calculate(mol)


def test_conformer_invalid_id(ref_mol):
    """
    Test the error raised of non-existing conformer id
    """
    calculator = spectrophore.SpectrophoreCalculator()

    with pytest.raises(ValueError, match="Conformer ID 2 is not a valid"):
        calculator.calculate(ref_mol, confID=2)


def test_conformer_2D(ref_mol):
    """
    Test the error raised of a 2D conformer
    """

    AllChem.Compute2DCoords(ref_mol)
    calculator = spectrophore.SpectrophoreCalculator()

    with pytest.raises(
        ValueError, match="The input molecule doesn't have a valid 3D conformation"
    ):
        calculator.calculate(ref_mol)


def test_spectrophore_run_1(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="all", stereo="none", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_1"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_2(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=10, normalization="all", stereo="none", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_2"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_3(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="none", stereo="none", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_3"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_4(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="mean", stereo="none", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_4"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_5(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="std", stereo="none", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_5"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_6(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="all", stereo="unique", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_6"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_7(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="all", stereo="mirror", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_7"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_8(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="all", stereo="all", resolution=3.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_8"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


def test_spectrophore_run_9(ref_mol, ref_data):
    """
    Test the output data
    """

    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=20, normalization="all", stereo="none", resolution=5.0
    )
    spec = calculator.calculate(ref_mol)

    ref_spec = ref_data["test_spectrophore_run_9"]

    np.testing.assert_allclose(spec, ref_spec, rtol=TOLERANCE)


@pytest.mark.parametrize("accuracy", [15, 20, 30])
@pytest.mark.parametrize("normalization", ["none", "mean", "std", "all"])
@pytest.mark.parametrize("stereo", ["none", "unique", "mirror", "all"])
@pytest.mark.parametrize("mode", ["classic", "full"])
def test_spectrophore_output(ref_mol, accuracy, normalization, stereo, mode):
    """
    Test the output dimensions and correct normalization with various spectrophore settings.
    """
    calculator = spectrophore.SpectrophoreCalculator(
        accuracy=accuracy, normalization=normalization, stereo=stereo, mode=mode
    )
    spec = calculator.calculate(ref_mol)

    if stereo == "none":
        num_probes = 12
    elif stereo == "unique":
        num_probes = 18
    elif stereo == "mirror":
        num_probes = 18
    elif stereo == "all":
        num_probes = 36

    # Check spec_size
    assert spec.shape[-1] == num_probes * 4

    # Check if number of frames are correct
    if mode == "full":
        assert spec.shape[0] == ((360 / accuracy) ** 2) * ((180 / accuracy) + 1)

    # Check normalization
    if normalization == "mean":
        spec = spec.reshape(-1, 4, num_probes)
        # Using atol when checking close to 0
        np.testing.assert_allclose(np.mean(spec, axis=-1), 0, atol=1e-5)

    elif normalization == "std":
        spec = spec.reshape(-1, 4, num_probes)
        np.testing.assert_allclose(np.std(spec, axis=-1), 1, rtol=TOLERANCE)

    elif normalization == "all":
        spec = spec.reshape(-1, 4, num_probes)
        # Using atol when checking close to 0
        np.testing.assert_allclose(np.mean(spec, axis=-1), 0, atol=1e-5)
        np.testing.assert_allclose(np.std(spec, axis=-1), 1, rtol=TOLERANCE)


def test_convert_sdf_alias(ref_data):
    """
    Test the convert sdf function using the 'spectrophore' alias.
    """

    input_path = os.path.join(os.path.dirname(__file__), "ref_mol.sdf")
    output_path = os.path.join(os.path.dirname(__file__), "output_test.csv")

    cmd_log = subprocess.run(
        f"spectrophore -i {input_path} -o {output_path} -p 1",
        shell=True,
        capture_output=True,
    )

    # Check successful output message
    assert "Successfully processed 1 out of 1 compounds" in cmd_log.stdout.decode(
        "utf-8"
    )

    # Check if output file is correct
    with open(output_path) as output_file:
        assert output_file.readline() == ref_data["test_convert_sdf"]

    # Remove file
    os.remove(output_path)


def test_convert_sdf_module(ref_data):
    """
    Test the convert sdf function using 'python -m spectrophore'.
    """

    input_path = os.path.join(os.path.dirname(__file__), "ref_mol.sdf")
    output_path = os.path.join(os.path.dirname(__file__), "output_test.csv")

    cmd_log = subprocess.run(
        f"python -m spectrophore -i {input_path} -o {output_path} -p 1",
        shell=True,
        capture_output=True,
    )

    # Check successful output message
    assert "Successfully processed 1 out of 1 compounds" in cmd_log.stdout.decode(
        "utf-8"
    )

    # Check if output file is correct
    with open(output_path) as output_file:
        assert output_file.readline() == ref_data["test_convert_sdf"]

    # Remove file
    os.remove(output_path)
