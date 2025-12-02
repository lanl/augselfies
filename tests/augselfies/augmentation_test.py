import pytest
import pandas as pd 
import rdkit.Chem 
import selfies as sf 
import augselfies.augmentation as aug 

@pytest.fixture
def df_dummy():
    return pd.DataFrame(data={"smiles":["C","O=C=O","Nc1ccccc1"]})


def test_equivalent_smiles():
    cyano_smile = "CN"
    equiv_smiles = aug.get_equivalent_smiles(cyano_smile,n_equivalents=2)
    assert "CN" in equiv_smiles, "original smiles string lost"
    assert "NC" in equiv_smiles, "secondary smiles not found"

def test_max_out_smiles_warning():
    cyano_smile = "ON"
    with pytest.warns():
        aug.get_equivalent_smiles(cyano_smile,n_equivalents=3), "Warning not thrown for too many equivalents"


def test_df_augment_size(df_dummy):
    n_equivs = 5
    len_init = len(df_dummy)
    aug_df = aug.augment_smiles_dataframe(df_dummy,n_equivalents=n_equivs,smiles_column='smiles')
    assert len_init*n_equivs == len(aug_df), "Augmented Dataframe wrong size"

def test_df_augment_selfies(df_dummy):
    aug_df = aug.augment_smiles_dataframe(df_dummy,n_equivalents=1,smiles_column='smiles',selfies_column='selfies')
    selfies_list = aug_df['selfies'].to_list()
    smiles_list = df_dummy['smiles'].to_list()
    canonical_conv_selfies = list(map(lambda x: rdkit.Chem.CanonSmiles(sf.decoder(x)),selfies_list))
    canonical_smiles = list(map(lambda x: rdkit.Chem.CanonSmiles(x),smiles_list))
    assert canonical_conv_selfies == canonical_smiles,"SELFIES conversion failed in augmentation"

if __name__ == "__main__":
    pytest.main()