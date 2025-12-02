import pytest
import rdkit.Chem
import selfies as sf 
import augselfies.numeralization

def update_mismatches(canon_smiles,canon_smiles_rt,mismatch_i,mismatch_o):
    if canon_smiles != canon_smiles_rt:
        mismatch_i.append(canon_smiles)
        mismatch_o.append(canon_smiles_rt)
def test_numselfies_roundtrip_basic(smiles_set):
    mismatch_i, mismatch_o = [], []
    for i, smiles in enumerate(smiles_set):
        canon_smiles = rdkit.Chem.CanonSmiles(smiles)
        num_selfies = augselfies.numeralization.selfies_to_num_selfies(sf.encoder(canon_smiles))
        canon_smiles_rt = rdkit.Chem.CanonSmiles(sf.decoder(augselfies.numeralization.num_selfies_to_selfies(num_selfies)))
        update_mismatches(canon_smiles,canon_smiles_rt,mismatch_i,mismatch_o)
    assert len(mismatch_i) == 0, f"Could not match {mismatch_i[0]}: got {mismatch_o[0]}"

def test_numselfies_roundtrip_basic_multiple_numeric_tokens(smiles_set):
    mismatch_i, mismatch_o = [], []
    for i, smiles in enumerate(smiles_set):
        canon_smiles = rdkit.Chem.CanonSmiles(smiles)
        num_selfies = augselfies.numeralization.selfies_to_num_selfies(sf.encoder(canon_smiles),single_numeric_token = False)
        canon_smiles_rt = rdkit.Chem.CanonSmiles(sf.decoder(augselfies.numeralization.num_selfies_to_selfies(num_selfies)))
        update_mismatches(canon_smiles,canon_smiles_rt,mismatch_i,mismatch_o)
    assert len(mismatch_i) == 0, f"Could not match {mismatch_i[0]}: got {mismatch_o[0]}"

def test_numselfies_roundtrip_basic_base_10(smiles_set):
    mismatch_i, mismatch_o = [], []
    for i, smiles in enumerate(smiles_set):
        canon_smiles = rdkit.Chem.CanonSmiles(smiles)
        num_selfies = augselfies.numeralization.selfies_to_num_selfies(sf.encoder(canon_smiles),base_int=10)
        canon_smiles_rt = rdkit.Chem.CanonSmiles(sf.decoder(augselfies.numeralization.num_selfies_to_selfies(num_selfies,base_int=10)))
        update_mismatches(canon_smiles,canon_smiles_rt,mismatch_i,mismatch_o)
    assert len(mismatch_i) == 0, f"Could not match {mismatch_i[0]}: got {mismatch_o[0]}"

