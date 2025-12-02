import warnings
from typing import Iterable, Optional, Union
import random
import pandas as pd 
from rdkit import Chem
import selfies as sf 

'''Routines for data augmentation of SMILES/SELFIES.'''

def get_equivalent_smiles(init_smiles:str,n_equivalents:int,trial_factor:int = 10)->list[str]:
    '''
    For an input SMILES string, returns a list of n equivalent strings.
    Args:
        init_smiles (str): 
        n_equivalents (int):
        trial_factor (int): Defaults to 10
            Multiplies the number of equivalents to attempt to create a new SMILES trial_factor*n_equivalents times.
            If n_equivalents cannot be made, method terminates with the already constructed amount. This typically happens for 
            small molecules for which n_equivalents may not exist. 
    Returns:
        list(str)
    '''
    total_tries = trial_factor*n_equivalents
    mol = Chem.MolFromSmiles(init_smiles)
    smiles_set = set()
    n_atoms = mol.GetNumAtoms()
    i = 0 
    while(n_equivalents > len(smiles_set)):
        i+=1
        if i > total_tries:
            warnings.warn(f"Could not produce {n_equivalents} SMILES after {total_tries}. Providing {len(smiles_set)} equivalents")
            break
        try:
            smiles_set.add(Chem.MolToSmiles(mol,canonical=False, rootedAtAtom=random.randint(0,n_atoms-1),doRandom=True))
        except:
            pass 
    return list(smiles_set)

def safe_selfies_encoder(smiles:str)->str:
    '''SMILES to SELFIES that passes empty string when encoding fails'''
    try:
        return sf.encoder(smiles)
    except:
        return ""

def safe_selfies_decoder(selfies:str)->str:
    '''SELFIES to SMILES that passes empty string when encoding fails'''
    try:
        return sf.decoder(selfies)
    except:
        return ""

def augment_smiles_dataframe(df:pd.DataFrame,n_equivalents:int, smiles_column:str = 'smiles',selfies_column:Optional[str]=None,
                             drop_null = True)->pd.DataFrame:
    '''
    Returns an augmented dataframe where the SMILES have been augmented by a factor of n_equivalents. 
    Note that in practice, the total number of strings after augmentation will be lower than the initial number of SMILES multiplied 
    by n_equivalents as n_equivalents cannot always be found. 
    For instance, augmenting a dataset of ~100k small molecules by n_equivalents = 5 will likely produce a total size of ~400k molecular strings
    after augmentation. 
    Args:
        df (pd.DataFrame): 
        n_equivalents (int): 
        smiles_column (str, optional): Defaults to 'smiles'.
        selfies_column (str, optional): Defaults to None
            If present, converts SMILES calculated to SELFIES.
    Returns:
        pd.DataFrame: augmented DataFrame
    '''
    aug_df = df.loc[df.index.repeat(n_equivalents)]
    aug_df.reset_index(inplace=True,drop = True)
    for i in range(len(df.index)):
        start_idx, end_idx = i*n_equivalents,(i+1)*n_equivalents-1
        cur_smiles = df[smiles_column].iloc[i]
        new_smiles_list = get_equivalent_smiles(cur_smiles,n_equivalents)
        if len(new_smiles_list)<n_equivalents:
            size_null = n_equivalents-len(new_smiles_list)
            new_smiles_list.extend([cur_smiles]*size_null)
        aug_df.loc[start_idx:end_idx,smiles_column] = new_smiles_list
    if selfies_column is not None:
        aug_df[selfies_column] = aug_df[smiles_column].map(lambda x: safe_selfies_encoder(x))
    if drop_null:
        aug_df.dropna()
    return aug_df