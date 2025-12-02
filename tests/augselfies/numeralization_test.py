import pytest
import rdkit.Chem 
import selfies 
from group_selfies import GroupGrammar
import group_selfies.constants
import augselfies.numeralization 


def test_branch_numeric():
    selfies_str = "[Branch1][Ring1]"
    num_selfies_str = augselfies.numeralization.selfies_to_num_selfies(selfies_str)
    num_tokens = augselfies.numeralization.get_token_list(num_selfies_str)
    assert num_tokens[1] == "[1]", "small branch failed"

def test_large_branch_numeric_hex():
    selfies_str = "[Branch2][Ring2][N]"
    num_selfies_str = augselfies.numeralization.selfies_to_num_selfies(selfies_str)
    num_tokens = augselfies.numeralization.get_token_list(num_selfies_str)
    assert num_tokens[1] == "[2a]", "large branch failed"

def test_reverse_large_branch_numeric():
    num_selfies_str = "[Branch][2a]"
    selfies_str = augselfies.numeralization.num_selfies_to_selfies(num_selfies_str)
    assert selfies_str == "[Branch2][Ring2][N]", "reverse large branch failed"

def test_ring_numeric():
    selfies_str = "[Ring1][O]"
    num_selfies_str = augselfies.numeralization.selfies_to_num_selfies(selfies_str)
    num_tokens = augselfies.numeralization.get_token_list(num_selfies_str)
    assert num_tokens[1] == "[9]", "small ring failed"

def test_aniline_rt():
    aniline_smiles = "Nc1ccccc1"
    aniline_selfies = selfies.encoder(aniline_smiles)
    aniline_num_selfies = augselfies.numeralization.selfies_to_num_selfies(aniline_selfies)
    aniline_selfies_rt = augselfies.numeralization.num_selfies_to_selfies(aniline_num_selfies)
    assert aniline_selfies == aniline_selfies_rt, "Aniline roundtrip failed"

def test_b12_rt():
    b12_smiles = r"CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)([O-])O[C@H](C)CNC(=O)CC[C@@]4([C@H]([C@@H]5[C@]6([C@@]([C@@H](/C(=C(/C7=N/C(=C\C8=N/C(=C(\C4=N5)/C)/[C@H](C8(C)C)CCC(=O)N)/[C@H]([C@]7(C)CC(=O)N)CCC(=O)N)\C)/[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3]"
    b12_selfies = selfies.encoder(b12_smiles)
    b12_num_selfies = augselfies.numeralization.selfies_to_num_selfies(b12_selfies)
    b12_selfies_rt = augselfies.numeralization.num_selfies_to_selfies(b12_num_selfies)
    assert b12_selfies == b12_selfies_rt , "B12 roundtrip failed"


def test_b12_smiles_rt():
    b12_smiles = r"CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)([O-])O[C@H](C)CNC(=O)CC[C@@]4([C@H]([C@@H]5[C@]6([C@@]([C@@H](/C(=C(/C7=N/C(=C\C8=N/C(=C(\C4=N5)/C)/[C@H](C8(C)C)CCC(=O)N)/[C@H]([C@]7(C)CC(=O)N)CCC(=O)N)\C)/[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3]"
    b12_selfies = selfies.encoder(b12_smiles)
    b12_smiles_rt = selfies.decoder(b12_selfies).replace("H1","H").replace("-1","-") #these are intended subsitutions
    assert b12_smiles == b12_smiles_rt, "B12 SMILES roundtrip failed"

def test_num_encoding():
    num_selfies_str = "[Ring][40][dummy]"
    new_tokens, num_values = augselfies.numeralization.num_selfies_to_text_num_encodings(num_selfies_str)
    assert new_tokens == ["[Ring]","[NUM]","[dummy]"]
    assert num_values == [1,64,1] #recall, this is a hex conversion by default

def test_num_encoding_decimal():
    num_selfies_str = "[Ring][32][dummy]"
    new_tokens, num_values = augselfies.numeralization.num_selfies_to_text_num_encodings(num_selfies_str,base_int=10)
    assert new_tokens == ["[Ring]","[NUM]","[dummy]"]
    assert num_values == [1,32,1] 

def test_num_encoding_omit_zero():
    num_selfies_str = "[Ring][0][dummy][Ring][2]"
    new_tokens, num_values = augselfies.numeralization.num_selfies_to_text_num_encodings(num_selfies_str,base_int = 16,omit_zero_token=True)
    assert new_tokens == ["[Ring]","[0]","[dummy]","[Ring]","[NUM]"]
    assert num_values == [1,1,1,1,2]  

def test_num_encoding_keep_zero():
    num_selfies_str = "[Ring][0][dummy][Ring][2]"
    new_tokens, num_values = augselfies.numeralization.num_selfies_to_text_num_encodings(num_selfies_str,base_int = 16,omit_zero_token=False)
    assert new_tokens == ["[Ring]","[NUM]","[dummy]","[Ring]","[NUM]"]
    assert num_values == [1,0,1,1,2]  

def test_num_selfies_large_ring_decimal():
    selfies = "[Ring2][Ring1][C]"
    num_selfies = augselfies.numeralization.selfies_to_num_selfies(selfies,base_int=10)
    assert num_selfies == "[Ring][16]", "num selfies for large ring failed "

def test_num_token_replacement():
    num_selfies = "[Ring1][e]"
    selfies = augselfies.numeralization.num_group_selfies_to_group_selfies(num_selfies,base_int=16)
    target_selfies = "[Ring1][S]"
    assert selfies == target_selfies, "num token replacement failed"

def test_num_token_replacement_base_10():
    num_selfies = "[:0dummy][16]"
    selfies = augselfies.numeralization.num_group_selfies_to_group_selfies(num_selfies,base_int=10)
    target_selfies = "[:0dummy][Ring1][C]"
    assert selfies == target_selfies, "num token replacement failed for base 10"


@pytest.fixture 
def float_to_convert():
    return 0

@pytest.fixture
def empty_to_convert():
    return ""

def test_num_selfies_to_selfies_float_handling(float_to_convert):
    selfies_float=  augselfies.numeralization.num_selfies_to_selfies(float_to_convert)
    assert selfies_float == "", "Float conversion did not produce empty string"

def test_selfies_to_num_selfies_float_handling(float_to_convert):
    selfies_float=  augselfies.numeralization.selfies_to_num_selfies(float_to_convert)
    assert selfies_float == "", "Float conversion did not produce empty string"


def test_num_selfies_to_selfies_empty_handling(empty_to_convert):
    selfies_float=  augselfies.numeralization.num_selfies_to_selfies(empty_to_convert)
    assert selfies_float == "", "Float conversion did not produce empty string"

def test_selfies_to_num_selfies_empty_handling(empty_to_convert):
    selfies_float=  augselfies.numeralization.selfies_to_num_selfies(empty_to_convert)
    assert selfies_float == "", "Float conversion did not produce empty string"

def test_hex_conversion():
    val = 1234
    hex_str = augselfies.numeralization.int_to_str(val,16)
    assert hex_str == hex(val)[2:], f"hex conversion of {val} failed"

def test_hex_conversion_negative():
    val = -1234
    hex_str = augselfies.numeralization.int_to_str(val,16)
    assert hex_str == "-"+hex(abs(val))[2:], f"hex conversion of {val} failed: Got {hex_str}"

def test_base_10_str_conversion():
    val = 623
    val_str = augselfies.numeralization.int_to_str(val,10)
    assert val_str == str(val), f"Base 10 str conversion of {val} failed: Got {val_str}"

def test_base_8_str_conversion():
    val = 10
    val_str = augselfies.numeralization.int_to_str(val,8)
    assert val_str == "12", f"Base 8 str conversion of {val} to 12 failed: Got{val_str}"

def test_base_8_str_conversion_neg():
    val = -10
    val_str = augselfies.numeralization.int_to_str(val,8)
    assert val_str == "-12", f"Base 8 str conversion of {val} to 12 failed: Got{val_str}"

def test_base_13_str_conversion():
    val = 25
    val_str = augselfies.numeralization.int_to_str(val,13)
    assert val_str == "1c", f"Base 13 conversion of 25 to 1c failed. Got {val_str}"

def test_merge_numeric_tokens():
    test_tokens = "[a][1][c][dummy][0]"
    return_tokens = augselfies.numeralization.merge_numeric_token_integers(test_tokens,base_int=16)
    assert return_tokens == "[a1c][dummy][0]"

def test_merge_numeric_tokens_base_10():
    test_tokens = "[5][3][1][dummy][0]"
    return_tokens = augselfies.numeralization.merge_numeric_token_integers(test_tokens,base_int=10)
    assert return_tokens == "[531][dummy][0]"

def test_merge_numeric_tokens_exception_toolarge():
    test_tokens = "[a][b][1][dummy][0]"
    with pytest.raises(ValueError):
        augselfies.numeralization.merge_numeric_token_integers(test_tokens,base_int=10)

if __name__ == "__main__":
    pytest.main()