import pytest
import augselfies.constants

def test_gselfies_default_index_size():
    assert len(augselfies.constants.DEFAULT_GROUPSELFIES_INDEX_ALPHABET) == 13

def test_selfies_default_index_size():
    assert len(augselfies.constants.DEFAULT_SELFIES_INDEX_ALPHABET) == 16

def test_num_selfies_default_index_size():
    assert len(augselfies.constants.DEFAULT_NUM_SELFIES_INDEX_ALPHABET) == 16