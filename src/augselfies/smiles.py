import multiprocessing.pool
from typing import Sequence
import multiprocessing
import rdkit.Chem
import selfies as sf 
import augselfies.numeralization

def get_canonical_smiles(chem_representation:None|str|rdkit.Chem.Mol)->str:
    '''
    Provides the canonical SMILES for a given representation (e.g., SMILES, SELFIES).
    Invalid transformations are returned as null string ""
    Args:
        chem_representation (None | str | rdkit.Chem.Mol):
    Returns:
        str: canonical SMILES or ""
    '''

    if chem_representation is None:
        mol = rdkit.Chem.MolFromSmiles("")
    if isinstance(chem_representation,str):
        smiles = ""
        try: 
            smiles = augselfies.numeralization.num_selfies_to_smiles(chem_representation)
        except:
            pass
        if smiles == "":
            try:
                smiles = sf.decoder(chem_representation)
            except:
                pass
        if smiles == "": 
            smiles = chem_representation
        try:
            mol = rdkit.Chem.MolFromSmiles(smiles)
        except:
            mol = rdkit.Chem.MolFromSmiles("")
    if isinstance(chem_representation,rdkit.Chem.Mol):
        mol = chem_representation
    try:
        canonical_smiles = rdkit.Chem.MolToSmiles(mol,canonical=True,kekuleSmiles = False)
    except:
        try:
            canonical_smiles = rdkit.Chem.MolToSmiles(mol,canonical=True)
        except:
            canonical_smiles = ""
    return canonical_smiles 

def get_canonical_smiles_list(chem_representations:Sequence[None|str|rdkit.Chem.Mol],pool: None|multiprocessing.pool.Pool = None)->list[str]:
    '''
    Provides the canonical SMILES for a given representation (e.g., SMILES, SELFIES).
    Invalid transformations are returned as null string ""
    Args:
        chem_representation (None | str | rdkit.Chem.Mol):
    Returns:
        list[str]: canonical SMILES or ""
    '''
    if pool is not None:
        smiles_list = list(pool.map(get_canonical_smiles,chem_representations))
    else:
        smiles_list = list(map(get_canonical_smiles,chem_representations))
    return smiles_list