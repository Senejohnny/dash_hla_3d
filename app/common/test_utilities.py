""" test file for utilities """
import pytest
from app.common.utilities import (
    get_hla_from_filename,
    flatten_dict_values
)


@pytest.mark.parametrize("test_input, expected_output",
    [
        ('DQA1_01_03-DQB1_06_01_V1.pdb', set(['DQA1*01:03', 'DQB1*06:01'])),
        ('A_33_01_V1.pdb', set(['A*33:01'])),
        ('DRB1_01_03_V1.pdb', set(['DRB1*01:03'])),
    ]
)
def test_get_hla_from_filename(test_input, expected_output): # pylint: disable=missing-function-docstring
    assert get_hla_from_filename(test_input) == expected_output

@pytest.mark.parametrize("test_input, expected_output",
    [
        ({'A': set('a'), 'B': set('b'), 'C': set('c')}, set(['a', 'b', 'c'])),
    ]
)
def test_flatten_dict_values(test_input, expected_output):
    assert flatten_dict_values(test_input) == expected_output
flatten_dict_values({'A': set('a'), 'B': set('b'), 'C': set('c')})
