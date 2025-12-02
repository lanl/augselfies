import pytest 
import rdkit.Chem
import selfies as sf
import augselfies.numeralization as num
import augselfies.smiles as smiles 
@pytest.fixture
def aniline_smiles():
    return "Nc1ccccc1"
@pytest.fixture
def aniline_selfies(aniline_smiles):
    return sf.encoder(aniline_smiles)

@pytest.fixture
def aniline_numselfies(aniline_selfies):
    return num.selfies_to_num_selfies(aniline_selfies)

@pytest.fixture
def aniline_mol(aniline_smiles):
    return rdkit.Chem.MolFromSmiles(aniline_smiles)

@pytest.fixture
def canonical_aniline_smiles(aniline_smiles):
    mol = rdkit.Chem.MolFromSmiles(aniline_smiles)
    return rdkit.Chem.MolToSmiles(mol,canonical=True)

def test_get_canonical_smiles_from_smiles(aniline_smiles,canonical_aniline_smiles):
    assert canonical_aniline_smiles == smiles.get_canonical_smiles(aniline_smiles), "Aniline SMILES did not canonicalize correctly"

def test_get_canonical_smiles_from_selfies(aniline_selfies,canonical_aniline_smiles):
    assert canonical_aniline_smiles == smiles.get_canonical_smiles(aniline_selfies), "Aniline SELFIES did not canonicalize correctly"

def test_get_canonical_smiles_from_numselfies(aniline_numselfies,canonical_aniline_smiles):
    assert canonical_aniline_smiles == smiles.get_canonical_smiles(aniline_numselfies), "Aniline numSELFIES did not canonicalize correctly"

def test_get_canonical_smiles_from_mol(aniline_mol,canonical_aniline_smiles):
    assert canonical_aniline_smiles == smiles.get_canonical_smiles(aniline_mol), "Aniline numSELFIES did not canonicalize correctly"

def test_get_canonical_smiles_from_null():
    assert "" == smiles.get_canonical_smiles(None), "None did not canonicalize correctly to ''"

def test_get_canonical_smiles_list(aniline_mol,canonical_aniline_smiles):
    aniline_list = [aniline_mol,aniline_mol]
    assert canonical_aniline_smiles == smiles.get_canonical_smiles_list(aniline_list)[1]