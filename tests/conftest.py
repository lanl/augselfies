import os
import pytest 
import pickle
import pandas as pd 
import rdkit.Chem 
from augselfies import PROJECT_ROOT
from group_selfies import Group 
from group_selfies.group_grammar import GroupGrammar


@pytest.fixture
def group_grammar_a_path():
    return os.path.join(PROJECT_ROOT,"tests","test_grammars","grammar_a")

@pytest.fixture
def group_grammar_a(group_grammar_a_path):
    with open(group_grammar_a_path,"rb") as grammar_file:
        group_grammar = pickle.load(grammar_file)
    return GroupGrammar(vocab=group_grammar.vocab)


@pytest.fixture
def chiral_grammar_path()->os.PathLike:
    return os.path.join(PROJECT_ROOT,"tests","test_grammars","chiral_grammar")

@pytest.fixture
def chiral_grammar(chiral_grammar_path)->GroupGrammar:
    with open(chiral_grammar_path,"rb") as grammar_file:
        group_grammar = pickle.load(grammar_file)
    return GroupGrammar(vocab=group_grammar.vocab)

@pytest.fixture
def smiles_set()->list[str]:
    df = pd.read_csv(os.path.join(PROJECT_ROOT,"tests","test_sets","smiles.csv"),delimiter = "\t")
    return df['smiles'].to_list()